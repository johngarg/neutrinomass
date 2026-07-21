#!/usr/bin/env python3

"""Bounded-memory exact deduplication for safe completion JSONL artifacts."""

from hashlib import sha256
from pathlib import Path
import resource
import sqlite3
import sys
import tempfile
from time import perf_counter

from neutrinomass.completions.completions import (
    are_equivalent_completions,
    exact_completion_bucket_key,
)
from neutrinomass.completions.equivalence import clear_interaction_graph_cache
from neutrinomass.completions.fingerprints import (
    completion_fingerprint,
    lagrangian_fingerprint,
)
from neutrinomass.database.serialization import (
    dumps_completion,
    iter_completion_jsonl,
    loads_completion,
)


def _file_digest(path):
    digest = sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _ordered_fingerprint_digest(connection, table):
    """Reproduce completion_fingerprint_digest using an on-disk sort."""

    digest = sha256()
    first = True
    for (fingerprint,) in connection.execute(
        f"SELECT fingerprint FROM {table} ORDER BY fingerprint"
    ):
        if not first:
            digest.update(b"\n")
        digest.update(fingerprint.encode("utf-8"))
        first = False
    return digest.hexdigest()


def _peak_memory_mib():
    maximum = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform != "darwin":
        maximum *= 1024
    return maximum / (1024 * 1024)


def _temporary_database_bytes(database_path):
    return sum(
        path.stat().st_size
        for path in (
            database_path,
            Path(str(database_path) + "-journal"),
            Path(str(database_path) + "-wal"),
        )
        if path.exists()
    )


def deduplicate_completion_jsonl(
    source,
    destination,
    *,
    work_dir=None,
    commit_interval=1000,
    representative_rank=None,
    maximum_rank=None,
):
    """Write exact classes from ``source`` to ``destination``.

    A SQLite index bounds the Python working set.  The Weisfeiler--Lehman
    Lagrangian fingerprint only prioritises comparisons within a physical
    bucket; every representative in that bucket remains eligible for exact
    contraction-graph isomorphism.  By default the first occurrence represents
    its class.  If ``representative_rank`` is supplied, an equivalent later
    occurrence with a larger integer rank replaces it.  Once ``maximum_rank``
    is reached, later equivalents need not be ranked.
    """

    if commit_interval < 1:
        raise ValueError("commit_interval must be positive")
    if maximum_rank is not None and representative_rank is None:
        raise ValueError("maximum_rank requires representative_rank")

    source = Path(source)
    destination = Path(destination)
    if source.resolve() == destination.resolve():
        raise ValueError("source and destination must differ")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if work_dir is not None:
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)

    started = perf_counter()
    input_records = 0
    retained_classes = 0
    candidate_hash_matches = 0
    exact_isomorphism_comparisons = 0
    representative_rank_evaluations = 0
    representative_replacements = 0
    maximum_database_bytes = 0

    temporary_output = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        prefix=f".{destination.name}.",
        suffix=".tmp",
        dir=destination.parent,
        delete=False,
    )
    temporary_output_path = Path(temporary_output.name)
    temporary_output.close()

    try:
        with tempfile.TemporaryDirectory(dir=work_dir) as directory:
            database_path = Path(directory) / "completion-deduplication.sqlite3"
            with sqlite3.connect(database_path) as connection:
                connection.execute("PRAGMA temp_store = FILE")
                connection.execute("PRAGMA cache_size = -8192")
                connection.execute("PRAGMA journal_mode = DELETE")
                connection.execute("PRAGMA synchronous = NORMAL")
                connection.executescript(
                    """
                    CREATE TABLE representatives (
                        sequence INTEGER PRIMARY KEY,
                        bucket_key TEXT NOT NULL,
                        candidate_hash TEXT NOT NULL,
                        payload TEXT NOT NULL,
                        fingerprint TEXT NOT NULL,
                        representative_rank INTEGER
                    );
                    CREATE INDEX representative_candidates
                        ON representatives(bucket_key, candidate_hash);
                    CREATE TABLE input_fingerprints (
                        fingerprint TEXT NOT NULL
                    );
                    """
                )

                clear_interaction_graph_cache()
                for candidate in iter_completion_jsonl(source):
                    try:
                        input_records += 1
                        bucket_key = repr(exact_completion_bucket_key(candidate))
                        candidate_hash = repr(lagrangian_fingerprint(candidate))
                        fingerprint = repr(completion_fingerprint(candidate))
                        connection.execute(
                            "INSERT INTO input_fingerprints(fingerprint) VALUES (?)",
                            (fingerprint,),
                        )

                        duplicate = False
                        possible_matches = connection.execute(
                            """
                            SELECT sequence, payload, candidate_hash,
                                   representative_rank
                            FROM representatives
                            WHERE bucket_key = ?
                            ORDER BY
                                CASE WHEN candidate_hash = ? THEN 0 ELSE 1 END,
                                sequence
                            """,
                            (bucket_key, candidate_hash),
                        )
                        for (
                            sequence,
                            payload,
                            known_hash,
                            known_rank,
                        ) in possible_matches:
                            if known_hash == candidate_hash:
                                candidate_hash_matches += 1
                            exact_isomorphism_comparisons += 1
                            representative = loads_completion(payload)
                            if are_equivalent_completions(
                                candidate, representative
                            ):
                                if representative_rank is not None:
                                    if known_rank is None:
                                        known_rank = representative_rank(
                                            representative
                                        )
                                        representative_rank_evaluations += 1
                                        connection.execute(
                                            """
                                            UPDATE representatives
                                            SET representative_rank = ?
                                            WHERE sequence = ?
                                            """,
                                            (known_rank, sequence),
                                        )
                                    if (
                                        maximum_rank is None
                                        or known_rank < maximum_rank
                                    ):
                                        candidate_rank = representative_rank(
                                            candidate
                                        )
                                        representative_rank_evaluations += 1
                                        if candidate_rank > known_rank:
                                            representative_replacements += 1
                                            connection.execute(
                                                """
                                                UPDATE representatives
                                                SET candidate_hash = ?,
                                                    payload = ?,
                                                    fingerprint = ?,
                                                    representative_rank = ?
                                                WHERE sequence = ?
                                                """,
                                                (
                                                    candidate_hash,
                                                    dumps_completion(candidate),
                                                    fingerprint,
                                                    candidate_rank,
                                                    sequence,
                                                ),
                                            )
                                duplicate = True
                                break

                        if not duplicate:
                            retained_classes += 1
                            connection.execute(
                                """
                                INSERT INTO representatives(
                                    sequence,
                                    bucket_key,
                                    candidate_hash,
                                    payload,
                                    fingerprint,
                                    representative_rank
                                ) VALUES (?, ?, ?, ?, ?, NULL)
                                """,
                                (
                                    input_records,
                                    bucket_key,
                                    candidate_hash,
                                    dumps_completion(candidate),
                                    fingerprint,
                                ),
                            )
                    finally:
                        clear_interaction_graph_cache()

                    if input_records % commit_interval == 0:
                        connection.commit()
                        maximum_database_bytes = max(
                            maximum_database_bytes,
                            _temporary_database_bytes(database_path),
                        )

                connection.commit()
                maximum_database_bytes = max(
                    maximum_database_bytes,
                    _temporary_database_bytes(database_path),
                )

                input_completion_digest = _ordered_fingerprint_digest(
                    connection, "input_fingerprints"
                )
                exact_completion_digest = _ordered_fingerprint_digest(
                    connection, "representatives"
                )
                coarse_buckets = connection.execute(
                    "SELECT COUNT(DISTINCT bucket_key) FROM representatives"
                ).fetchone()[0]
                maximum_bucket_size = connection.execute(
                    """
                    SELECT COALESCE(MAX(bucket_size), 0)
                    FROM (
                        SELECT COUNT(*) AS bucket_size
                        FROM representatives
                        GROUP BY bucket_key
                    )
                    """
                ).fetchone()[0]

                with temporary_output_path.open("w", encoding="utf-8") as stream:
                    for (payload,) in connection.execute(
                        "SELECT payload FROM representatives ORDER BY sequence"
                    ):
                        stream.write(payload + "\n")

            temporary_output_path.replace(destination)

        return {
            "source": str(source.resolve()),
            "destination": str(destination.resolve()),
            "input_records": input_records,
            "exact_classes": retained_classes,
            "coarse_buckets": coarse_buckets,
            "maximum_bucket_size": maximum_bucket_size,
            "candidate_hash_matches": candidate_hash_matches,
            "exact_isomorphism_comparisons": exact_isomorphism_comparisons,
            "representative_rank_evaluations": (
                representative_rank_evaluations
            ),
            "representative_replacements": representative_replacements,
            "input_completion_digest": input_completion_digest,
            "exact_completion_digest": exact_completion_digest,
            "source_sha256": _file_digest(source),
            "destination_sha256": _file_digest(destination),
            "temporary_database_peak_bytes": maximum_database_bytes,
            "process_peak_memory_mib": _peak_memory_mib(),
            "wall_time_seconds": perf_counter() - started,
        }
    finally:
        temporary_output_path.unlink(missing_ok=True)
