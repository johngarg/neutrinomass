"""Exact reduction that checkpoints inside a large equivalence bucket."""

import gzip
import json
from pathlib import Path

from neutrinomass.completions.completions import (
    are_equivalent_completions, exact_completion_bucket_key,
)
from neutrinomass.completions.fingerprints import completion_fingerprint
from neutrinomass.database.checkpoints import (
    CheckpointStore, WorkBudget, digest_file,
)
from neutrinomass.database.deduplication import _ordered_fingerprint_digest
from neutrinomass.database.global_equivalence import DecodedCache
from neutrinomass.database.serialization import loads_completion


def lines_from_paths(paths):
    for path in paths:
        path = Path(path)
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path,"rt",encoding="utf-8") as stream:
            for line in stream:
                if line.strip():
                    yield line.rstrip("\n")


def input_contract(paths):
    return [dict(sha256=digest_file(path),bytes=Path(path).stat().st_size) for path in paths]


def reduce_resumable(paths, directory, contract, *, budget=None, scratch=None,
                     representative_rank=None, maximum_rank=None, relation="physical"):
    """Preserve first-occurrence order, exact matching and survivor ranking.

    A candidate hash only orders work. The saved comparison cursor contains
    both the hash-priority group and representative sequence, so a restart
    cannot omit the remainder of a bucket or repeat its completed prefix.
    """
    paths = tuple(Path(p) for p in paths)
    budget = budget or WorkBudget()
    if relation == "physical":
        bucket_key, equivalent = exact_completion_bucket_key, are_equivalent_completions
    elif relation == "global-lagrangian":
        from neutrinomass.database.global_equivalence import global_lagrangian_bucket, equivalent_global_lagrangians
        bucket_key, equivalent = global_lagrangian_bucket, equivalent_global_lagrangians
    else:
        raise ValueError("Unknown exact-equivalence relation")
    identity = dict(contract, algorithm="exact-reduction-v1", relation=relation, inputs=input_contract(paths),
        ranking=None if representative_rank is None else representative_rank.__name__,
        maximum_rank=maximum_rank)
    if maximum_rank is not None and representative_rank is None:
        raise ValueError("maximum_rank requires representative_rank")
    with CheckpointStore(directory,identity) as store:
        if store.state["complete"]:
            return store.state
        with store.database(scratch) as db:
            db.executescript("""
                CREATE TABLE IF NOT EXISTS representatives (
                    sequence INTEGER PRIMARY KEY, bucket_key TEXT NOT NULL,
                    candidate_hash TEXT NOT NULL, payload TEXT NOT NULL,
                    fingerprint TEXT NOT NULL, representative_rank INTEGER);
                CREATE INDEX IF NOT EXISTS representative_buckets
                    ON representatives(bucket_key,candidate_hash);
                CREATE TABLE IF NOT EXISTS input_fingerprints (
                    sequence INTEGER PRIMARY KEY, fingerprint TEXT NOT NULL);
            """)
            progress = dict(processed=0, comparisons=0, replacements=0, pending=None,
                            **store.state["progress"]) if not store.state["progress"] else dict(store.state["progress"])
            cache = DecodedCache()

            def save():
                store.snapshot(db,progress)

            for ordinal, payload in enumerate(lines_from_paths(paths),1):
                if ordinal <= progress["processed"]:
                    continue
                candidate = loads_completion(payload)
                pending = progress["pending"]
                if pending is None:
                    fingerprint = completion_fingerprint(candidate)
                    candidate_hash = repr(fingerprint[4])
                    bucket = repr(bucket_key(candidate))
                    pending = dict(ordinal=ordinal, payload=payload, bucket=bucket,
                        candidate_hash=candidate_hash, fingerprint=repr(fingerprint), cursor=(-1,-1))
                    progress["pending"] = pending
                    db.execute("INSERT INTO input_fingerprints VALUES (?,?)",(ordinal,repr(fingerprint)))
                elif pending["ordinal"] != ordinal or pending["payload"] != payload:
                    raise ValueError("Input sequence changed during reduction resume")
                duplicate = False
                priority, sequence = pending["cursor"]
                rows = db.execute("""
                        SELECT sequence,payload,representative_rank,
                               CASE WHEN candidate_hash=? THEN 0 ELSE 1 END AS priority
                        FROM representatives WHERE bucket_key=?
                        AND (CASE WHEN candidate_hash=? THEN 0 ELSE 1 END > ?
                             OR (CASE WHEN candidate_hash=? THEN 0 ELSE 1 END = ? AND sequence > ?))
                        ORDER BY priority,sequence
                    """,(pending["candidate_hash"],pending["bucket"],pending["candidate_hash"],priority,
                          pending["candidate_hash"],priority,sequence))
                for row in rows:
                    sequence, known_payload, known_rank, priority = row
                    representative = cache.get(known_payload)
                    matched = equivalent(candidate,representative)
                    progress["comparisons"] += 1
                    pending["cursor"] = (priority,sequence)
                    if matched:
                        if representative_rank is not None:
                            if known_rank is None:
                                known_rank = representative_rank(representative)
                                db.execute("UPDATE representatives SET representative_rank=? WHERE sequence=?",
                                           (known_rank,sequence))
                            if maximum_rank is None or known_rank < maximum_rank:
                                rank = representative_rank(candidate)
                                if rank > known_rank:
                                    db.execute("""UPDATE representatives SET candidate_hash=?,payload=?,
                                        fingerprint=?,representative_rank=? WHERE sequence=?""",
                                        (pending["candidate_hash"],payload,pending["fingerprint"],rank,sequence))
                                    progress["replacements"] += 1
                        duplicate = True
                        # Finish this input atomically before checkpointing an
                        # accepted match; otherwise a resume would insert it.
                        break
                    budget.boundary(save,progress=dict(stage="exact-reduction",
                        processed=progress["processed"],comparisons=progress["comparisons"]))
                if not duplicate:
                    db.execute("INSERT INTO representatives VALUES (?,?,?,?,?,NULL)",
                        (ordinal,pending["bucket"],pending["candidate_hash"],payload,pending["fingerprint"]))
                progress["processed"] = ordinal
                progress["pending"] = None
                budget.boundary(save,progress=dict(stage="exact-reduction",
                    processed=ordinal,comparisons=progress["comparisons"]))
            if progress["pending"] is not None:
                raise ValueError("Input ended before the pending reduction candidate")
            output = store.publish_lines(payload for (payload,) in db.execute(
                "SELECT payload FROM representatives ORDER BY sequence"))
            progress.update(
                input_records=progress["processed"], exact_classes=output["records"],
                input_completion_digest=_ordered_fingerprint_digest(db,"input_fingerprints"),
                exact_completion_digest=_ordered_fingerprint_digest(db,"representatives"),
                coarse_buckets=db.execute("SELECT COUNT(DISTINCT bucket_key) FROM representatives").fetchone()[0])
            store.save(progress,artifacts={"output":output})
            store.snapshot(db,progress,complete=True)
            return store.state
