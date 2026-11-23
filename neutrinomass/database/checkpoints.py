"""Durable, checksum-verified stage checkpoints and cooperative time budgets.

Mutable SQLite work is copied from a sealed snapshot into local scratch.
Publication orders immutable outputs before the atomic manifest. An abrupt
exit can leave unreferenced files, but cannot certify a partial output.
"""

from contextlib import contextmanager
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import sqlite3
import sys
import tempfile
import time
import threading


SCHEMA_VERSION = 1


def digest_file(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def identity_hash(value):
    return hashlib.sha256(canonical_json(value).encode()).hexdigest()


def sync_directory(path):
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".publish-", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as stream:
            stream.write(canonical_json(value)+"\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        sync_directory(path.parent)
    finally:
        Path(temporary).unlink(missing_ok=True)


def memory_sample():
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result = {"peak_rss_mib": peak/(1024**2 if sys.platform == "darwin" else 1024)}
    try:
        for line in Path("/proc/self/status").read_text().splitlines():
            if line.startswith("VmRSS:"):
                result["rss_mib"] = int(line.split()[1])/1024
        relative = next(line.split(":",2)[2] for line in Path("/proc/self/cgroup").read_text().splitlines()
                        if line.startswith("0::"))
        root = Path("/sys/fs/cgroup")/relative.lstrip("/")
        result["cgroup_memory"] = {line.split()[0]: int(line.split()[1])
            for line in (root/"memory.stat").read_text().splitlines()
            if line.split()[0] in {"anon","file","shmem"}}
        result["cgroup_events"] = {k:int(v) for k,v in
            (line.split() for line in (root/"memory.events").read_text().splitlines())}
    except (OSError, ValueError, StopIteration):
        pass
    return result


class BudgetExpired(Exception):
    """Exit 75 means a saved continuation, not a completed scientific task."""


class WorkBudget:
    def __init__(self, seconds=None, *, checkpoint_seconds=300, fault_after=None):
        self.started = time.monotonic()
        self.deadline = self.started + seconds if seconds is not None else float("inf")
        self.checkpoint_seconds = checkpoint_seconds
        self.last_checkpoint = self.started
        self.last_progress = 0
        self.stop_requested = False
        self.fault_after = fault_after
        self.work = 0
        self.last_details = {}
        self.scratch_root = None

    def scratch_sample(self):
        if self.scratch_root is None:
            return {}
        total=0
        for path in Path(self.scratch_root).rglob("*"):
            try:
                if path.is_file():
                    total+=path.stat().st_size
            except FileNotFoundError:
                pass
        return {"scratch_bytes":total}

    @contextmanager
    def signals(self):
        previous = {}
        finished = threading.Event()
        def heartbeat():
            while not finished.wait(60):
                print(canonical_json(dict(event="worker_heartbeat",
                    elapsed_seconds=time.monotonic()-self.started,
                    checkpoint_age_seconds=time.monotonic()-self.last_checkpoint,
                    stop_requested=self.stopping, work=self.work,
                    **memory_sample(), **self.scratch_sample(), **self.last_details)),flush=True)
        thread = threading.Thread(target=heartbeat,daemon=True)
        def request_stop(signum, frame):
            self.stop_requested = True
        for sig in (signal.SIGUSR1, signal.SIGTERM):
            previous[sig] = signal.signal(sig, request_stop)
        thread.start()
        try:
            yield self
        finally:
            finished.set()
            thread.join(timeout=1)
            for sig, handler in previous.items():
                signal.signal(sig, handler)

    @property
    def stopping(self):
        return self.stop_requested or time.monotonic() >= self.deadline

    def boundary(self, save, *, progress=None, force=False):
        self.work += 1
        self.last_details = dict(progress or {})
        now = time.monotonic()
        if self.fault_after is not None and self.work >= self.fault_after:
            self.stop_requested = True
        if now-self.last_progress >= 60:
            print(canonical_json(dict(event="stage_progress", elapsed_seconds=now-self.started,
                work=self.work, **memory_sample(), **(progress or {}))), flush=True)
            self.last_progress = now
        if force or self.stopping or now-self.last_checkpoint >= self.checkpoint_seconds:
            save()
            self.last_checkpoint = time.monotonic()
        if self.stopping:
            raise BudgetExpired("Checkpoint saved; time budget or signal requested continuation")


class CheckpointStore:
    def __init__(self, directory, contract):
        self.directory = Path(directory)
        self.contract = json.loads(canonical_json(contract))
        self.state = None
        self.lock = None

    def __enter__(self):
        self.last_saved=time.monotonic()
        self.directory.mkdir(parents=True, exist_ok=True)
        self.lock = (self.directory/"worker.lock").open("a+")
        try:
            fcntl.flock(self.lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            manifest = self.directory/"checkpoint.json"
            if manifest.exists():
                state = json.loads(manifest.read_text())
                if state.get("schema_version") != SCHEMA_VERSION or state.get("contract") != self.contract:
                    raise ValueError("Incompatible checkpoint contract")
                if state.get("contract_sha256") != identity_hash(self.contract):
                    raise ValueError("Damaged checkpoint identity")
                for artifact in state.get("artifacts", {}).values():
                    self.verify(artifact)
                self.state = state
            else:
                self.state = dict(schema_version=SCHEMA_VERSION, contract=self.contract,
                    contract_sha256=identity_hash(self.contract), complete=False, artifacts={}, progress={})
            return self
        except BaseException:
            self.lock.close()
            raise

    def __exit__(self, *args):
        fcntl.flock(self.lock, fcntl.LOCK_UN)
        self.lock.close()

    def verify(self, artifact):
        path = self.directory/artifact["path"]
        if path.parent.resolve() != self.directory.resolve():
            raise ValueError("Checkpoint artifact escapes its stage directory")
        if path.stat().st_size != artifact["bytes"] or digest_file(path) != artifact["sha256"]:
            raise ValueError(f"Damaged checkpoint artifact: {path}")
        return path

    def publish_file(self, source, *, suffix):
        source = Path(source)
        digest = digest_file(source)
        destination = self.directory/(digest+suffix)
        if destination.exists():
            if digest_file(destination) != digest:
                raise ValueError("Existing immutable artifact is corrupt")
        else:
            # source can reside on another filesystem, so copy before rename.
            import shutil
            fd, name = tempfile.mkstemp(prefix=".artifact-", dir=self.directory)
            try:
                with source.open("rb") as src, os.fdopen(fd,"wb") as dst:
                    shutil.copyfileobj(src, dst, 1024**2)
                    dst.flush(); os.fsync(dst.fileno())
                if digest_file(name) != digest:
                    raise ValueError("Artifact changed while publishing")
                os.replace(name, destination)
                sync_directory(self.directory)
            finally:
                Path(name).unlink(missing_ok=True)
        return dict(path=destination.name, sha256=digest, bytes=destination.stat().st_size)

    def publish_lines(self, lines):
        fd, name = tempfile.mkstemp(prefix=".chunk-", dir=self.directory)
        count = 0
        plain = hashlib.sha256()
        try:
            with os.fdopen(fd,"wb") as dest:
                with gzip.GzipFile(fileobj=dest, mode="wb", filename="", mtime=0, compresslevel=2) as zipped:
                    for line in lines:
                        data = (line.rstrip("\n")+"\n").encode()
                        zipped.write(data); plain.update(data); count += 1
                dest.flush(); os.fsync(dest.fileno())
            with gzip.open(name,"rb") as restored:
                if hashlib.file_digest(restored,"sha256").hexdigest() != plain.hexdigest():
                    raise ValueError("Compressed chunk round-trip changed")
            artifact = self.publish_file(name, suffix=".jsonl.gz")
            artifact.update(records=count, uncompressed_sha256=plain.hexdigest())
            return artifact
        finally:
            Path(name).unlink(missing_ok=True)

    def save(self, progress, *, artifacts=None, complete=False):
        retained = dict(self.state["artifacts"])
        retained.update(artifacts or {})
        self.state = dict(self.state, progress=progress, artifacts=retained,
            complete=complete, updated_unix=time.time(), memory=memory_sample(),
            active_seconds=self.state.get("active_seconds",0)+time.monotonic()-self.last_saved)
        atomic_json(self.directory/"checkpoint.json", self.state)
        self.last_saved=time.monotonic()

    def snapshot(self, connection, progress, *, complete=False):
        connection.commit()
        fd, name = tempfile.mkstemp(prefix=".snapshot-", dir=self.directory)
        os.close(fd)
        try:
            with sqlite3.connect(name) as copy:
                connection.backup(copy)
                if copy.execute("PRAGMA quick_check").fetchone()[0] != "ok":
                    raise ValueError("Invalid SQLite snapshot")
            with open(name,"rb") as stream:
                os.fsync(stream.fileno())
            old = self.state["artifacts"].get("database")
            artifact = self.publish_file(name, suffix=".sqlite3")
            self.save(progress, artifacts={"database": artifact}, complete=complete)
            # Only the previous sealed database is superseded. Raw/output
            # chunks and every artifact outside this stage are untouched.
            if old and old["path"] != artifact["path"]:
                (self.directory/old["path"]).unlink(missing_ok=True)
        finally:
            Path(name).unlink(missing_ok=True)

    @contextmanager
    def database(self, scratch=None):
        with tempfile.TemporaryDirectory(prefix="census-", dir=scratch) as directory:
            path = Path(directory)/"work.sqlite3"
            if "database" in self.state["artifacts"]:
                import shutil
                shutil.copyfile(self.verify(self.state["artifacts"]["database"]),path)
            connection = sqlite3.connect(path)
            connection.execute("PRAGMA temp_store=FILE")
            connection.execute("PRAGMA cache_size=-8192")
            try:
                yield connection
            finally:
                connection.close()
