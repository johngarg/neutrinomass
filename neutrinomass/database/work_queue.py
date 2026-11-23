"""Persistent dependency queue; computation remains in Slurm workers."""

from contextlib import contextmanager
import json
from pathlib import Path
import sqlite3
import time

from neutrinomass.database.checkpoints import canonical_json,identity_hash


class WorkQueue:
    def __init__(self,path):
        self.path=Path(path)
        self.path.parent.mkdir(parents=True,exist_ok=True)
        with self.connection() as db:
            db.executescript("""
                CREATE TABLE IF NOT EXISTS units (
                    id TEXT PRIMARY KEY,phase TEXT,seed INTEGER,payload TEXT,status TEXT,
                    worker TEXT,progress TEXT,no_progress INTEGER DEFAULT 0,
                    infrastructure_retries INTEGER DEFAULT 0,result TEXT,diagnostic TEXT,
                    updated REAL);
                CREATE TABLE IF NOT EXISTS dependencies (unit TEXT,dependency TEXT,
                    PRIMARY KEY(unit,dependency));
                CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY,value TEXT);
            """)

    @contextmanager
    def connection(self):
        db=sqlite3.connect(self.path,timeout=60)
        db.row_factory=sqlite3.Row
        try:
            yield db
            db.commit()
        except BaseException:
            db.rollback()
            raise
        finally:
            db.close()

    def add(self,identifier,payload,*,phase,seed,dependencies=()):
        with self.connection() as db:
            db.execute("BEGIN IMMEDIATE")
            old=db.execute("SELECT payload,phase,seed FROM units WHERE id=?",(identifier,)).fetchone()
            encoded=canonical_json(payload)
            if old is not None:
                if tuple(old) != (encoded,phase,seed):
                    raise ValueError("Stable work identifier has incompatible input")
                actual={row[0] for row in db.execute("SELECT dependency FROM dependencies WHERE unit=?",(identifier,))}
                if actual != set(dependencies):
                    raise ValueError("Stable work dependencies changed")
                return
            db.execute("INSERT INTO units(id,phase,seed,payload,status,updated) VALUES (?,?,?,?,?,?)",
                       (identifier,phase,seed,encoded,"pending",time.time()))
            db.executemany("INSERT INTO dependencies VALUES (?,?)",((identifier,d) for d in dependencies))

    def claim(self,worker,*,phase,seed):
        with self.connection() as db:
            db.execute("BEGIN IMMEDIATE")
            if db.execute("SELECT 1 FROM units WHERE status='running' AND worker=?",(worker,)).fetchone():
                raise ValueError("Worker already owns an unfinished unit")
            row=db.execute("""SELECT * FROM units u WHERE status='pending' AND phase=? AND seed=?
                AND NOT EXISTS (SELECT 1 FROM dependencies d LEFT JOIN units p ON p.id=d.dependency
                    WHERE d.unit=u.id AND (p.status IS NULL OR p.status!='complete'))
                ORDER BY u.rowid LIMIT 1""",(phase,seed)).fetchone()
            if row is None:
                return None
            db.execute("UPDATE units SET status='running',worker=?,updated=? WHERE id=?",(worker,time.time(),row["id"]))
            return dict(id=row["id"],payload=json.loads(row["payload"]),phase=phase,seed=seed)

    def finish(self,identifier,worker,result):
        encoded=canonical_json(result)
        with self.connection() as db:
            db.execute("BEGIN IMMEDIATE")
            row=db.execute("SELECT * FROM units WHERE id=?",(identifier,)).fetchone()
            if row["status"] == "complete" and row["result"] == encoded:
                return
            if row["status"] != "running" or row["worker"] != worker:
                raise ValueError("Only the owning worker may publish a unit")
            db.execute("UPDATE units SET status='complete',result=?,updated=? WHERE id=?",(encoded,time.time(),identifier))

    def continuation(self,identifier,worker,progress):
        signature=identity_hash(progress)
        with self.connection() as db:
            db.execute("BEGIN IMMEDIATE")
            row=db.execute("SELECT * FROM units WHERE id=?",(identifier,)).fetchone()
            if row["status"] != "running" or row["worker"] != worker:
                raise ValueError("Continuation is not owned by this worker")
            repetitions=row["no_progress"]+1 if row["progress"] == signature else 0
            status="blocked" if repetitions >= 2 else "pending"
            diagnostic="REPEATED_NO_PROGRESS" if status == "blocked" else None
            db.execute("UPDATE units SET status=?,progress=?,no_progress=?,diagnostic=?,updated=? WHERE id=?",
                       (status,signature,repetitions,diagnostic,time.time(),identifier))
            return status

    def fail(self,identifier,worker,diagnostic):
        with self.connection() as db:
            changed=db.execute("UPDATE units SET status='blocked',diagnostic=?,updated=? WHERE id=? AND worker=? AND status='running'",
                       (canonical_json(diagnostic),time.time(),identifier,worker)).rowcount
            if changed != 1:
                raise ValueError("Failure is not owned by this worker")

    def recover_infrastructure(self,worker,reason):
        if reason not in {"NODE_FAIL","PREEMPTED","BOOT_FAIL"}:
            raise ValueError("Failure is not an identified infrastructure condition")
        with self.connection() as db:
            db.execute("BEGIN IMMEDIATE")
            for row in db.execute("SELECT id,infrastructure_retries FROM units WHERE worker=? AND status='running'",(worker,)).fetchall():
                retries=row["infrastructure_retries"]+1
                status="pending" if retries <= 2 else "blocked"
                db.execute("UPDATE units SET status=?,infrastructure_retries=?,diagnostic=?,updated=? WHERE id=?",
                    (status,retries,reason,time.time(),row["id"]))

    def summary(self):
        with self.connection() as db:
            return dict(counts=[dict(row) for row in db.execute(
                "SELECT phase,seed,status,COUNT(*) AS units FROM units GROUP BY phase,seed,status")],
                blocked=[dict(row) for row in db.execute("SELECT id,diagnostic FROM units WHERE status='blocked'")])


def progress_signature(directory):
    root=Path(directory)
    return [(str(path.relative_to(root)),data["progress"],data["complete"])
            for path in sorted(root.rglob("checkpoint.json"))
            for data in [json.loads(path.read_text())]]
