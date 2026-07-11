#!/usr/bin/env python3

"""Reproducible, bounded-memory building blocks for database regeneration."""

from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations
import json
import math
import os
from pathlib import Path
import resource
import sqlite3
import sys
import tempfile
from time import perf_counter

from neutrinomass.completions.completions import (
    HistoricalDerivativeBasis,
    PROJECTED_LORENTZ_OPERATORS,
    append_unique_completions,
    are_equivalent_completions,
    base_exotic_label,
    deriv_operator_completions,
    exact_completion_bucket_key,
    is_singlet,
    operator_completions,
    operator_strip_derivs,
)
from neutrinomass.completions.equivalence import clear_interaction_graph_cache
from neutrinomass.completions.fingerprints import (
    completion_fingerprint,
    democratic_model_fingerprint,
)
from neutrinomass.completions.operators import (
    DERIV_EFF_OPERATORS,
    EFF_OPERATORS,
)
from neutrinomass.database.closures import (
    neutrino_mass_estimate,
    numerical_np_scale_estimate,
)
from neutrinomass.database.database import read_completions
from neutrinomass.database.deduplication import deduplicate_completion_jsonl
from neutrinomass.database.serialization import (
    dumps_completion,
    iter_completion_jsonl,
    loads_completion,
)


REGULAR_KIND = "regular"
DERIVATIVE_KIND = "derivative"


def file_sha256(path):
    digest = sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def peak_memory_mib():
    maximum = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform != "darwin":
        maximum *= 1024
    return maximum / (1024 * 1024)


def operator_registry():
    overlap = set(EFF_OPERATORS) & set(DERIV_EFF_OPERATORS)
    if overlap:
        raise ValueError(f"Operator registries overlap: {sorted(overlap)}")
    return {**EFF_OPERATORS, **DERIV_EFF_OPERATORS}


def operator_kind(operator_name):
    if operator_name in EFF_OPERATORS:
        return REGULAR_KIND
    if operator_name in DERIV_EFF_OPERATORS:
        return DERIVATIVE_KIND
    raise KeyError(operator_name)


def operator_inventory(legacy_dir):
    """Return the exact registry-to-legacy-file mapping or fail visibly."""

    legacy_dir = Path(legacy_dir)
    legacy = {
        path.stem.removeprefix("op_"): path
        for path in legacy_dir.glob("op_*.dat")
    }
    registry = operator_registry()
    missing_definitions = sorted(set(legacy) - set(registry))
    missing_legacy = sorted(set(registry) - set(legacy))
    if missing_definitions or missing_legacy:
        raise ValueError(
            "Registry/legacy inventory mismatch: "
            f"missing definitions={missing_definitions}, "
            f"missing legacy files={missing_legacy}"
        )
    return tuple(
        {
            "operator": name,
            "kind": operator_kind(name),
            "derivatives": (
                operator_strip_derivs(registry[name].operator)["n_derivs"]
                if name in DERIV_EFF_OPERATORS
                else 0
            ),
            "legacy_path": str(legacy[name].resolve()),
            "legacy_sha256": file_sha256(legacy[name]),
        }
        for name in sorted(registry)
    )


def completion_stream(operator_name):
    if operator_name in EFF_OPERATORS:
        return operator_completions(EFF_OPERATORS[operator_name])
    if operator_name in DERIV_EFF_OPERATORS:
        return iter(deriv_operator_completions(DERIV_EFF_OPERATORS[operator_name]))
    raise KeyError(operator_name)


def topology_key(completion):
    return f"{completion.topology}|{completion.canonical_topology}"


def model_strings(completion):
    def quantum_number_string(info):
        lorentz, colour_up, colour_down, isospin, (_, baryon), (_, hypercharge) = info
        return (
            f"{lorentz},{colour_up}{colour_down},{isospin},"
            f"{hypercharge},{baryon}"
        )

    return tuple(
        sorted(
            quantum_number_string(info)
            for info in democratic_model_fingerprint(completion)
        )
    )


def validate_completion(completion):
    if any(term.safe_simplify() == 0 for term in completion.terms):
        raise ValueError("vanishing UV interaction")
    if any(
        not is_singlet(term)
        or sum(field.mass_dim for field in term.fields) > 4
        for term in completion.terms
    ):
        raise ValueError("non-renormalisable or non-singlet UV interaction")

    if not completion.derivative_routes:
        return
    if len(completion.derivative_routes) != 1:
        raise ValueError("routed completion does not have exactly one route")
    route = completion.derivative_routes[0]
    if not completion.graph.has_edge(*route.edge):
        raise ValueError(f"recorded route edge is absent: {route.edge}")
    particle = completion.graph.edges[route.edge]["particle"]
    if base_exotic_label(particle) != base_exotic_label(route.numerator_field):
        raise ValueError("route numerator does not match the edge particle")


def round_trip_signature(completion):
    projection = getattr(completion, "lorentz_projection", None)
    return (
        completion_fingerprint(completion),
        tuple(completion.derivative_routes),
        tuple(projection) if projection is not None else None,
    )


class OrderedSignatureDigest:
    def __init__(self):
        self._digest = sha256()
        self._first = True

    def update(self, value):
        if not self._first:
            self._digest.update(b"\n")
        self._digest.update(repr(value).encode("utf-8"))
        self._first = False

    def hexdigest(self):
        return self._digest.hexdigest()


def _new_stats():
    return {
        "records": 0,
        "local": 0,
        "routed": 0,
        "topologies": Counter(),
    }


def _update_stats(stats, completion):
    stats["records"] += 1
    routed = bool(completion.derivative_routes)
    stats["routed" if routed else "local"] += 1
    stats["topologies"][topology_key(completion)] += 1


def _serialise_stats(stats):
    return {
        "records": stats["records"],
        "local": stats["local"],
        "routed": stats["routed"],
        "topologies": dict(sorted(stats["topologies"].items())),
    }


def write_generated_artifact(path, completions):
    """Atomically write, validate and stream-audit generated completions."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    source_stats = _new_stats()
    source_digest = OrderedSignatureDigest()
    temporary = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
        delete=False,
    )
    temporary_path = Path(temporary.name)
    try:
        with temporary:
            for completion in completions:
                validate_completion(completion)
                _update_stats(source_stats, completion)
                source_digest.update(round_trip_signature(completion))
                temporary.write(dumps_completion(completion) + "\n")
        temporary_path.replace(path)
    finally:
        temporary_path.unlink(missing_ok=True)

    decoded_stats = _new_stats()
    decoded_digest = OrderedSignatureDigest()
    for completion in iter_completion_jsonl(path):
        validate_completion(completion)
        _update_stats(decoded_stats, completion)
        decoded_digest.update(round_trip_signature(completion))

    source = _serialise_stats(source_stats)
    decoded = _serialise_stats(decoded_stats)
    if source != decoded:
        raise ValueError("record statistics changed across JSONL round trip")
    if source_digest.hexdigest() != decoded_digest.hexdigest():
        raise ValueError("completion metadata changed across JSONL round trip")
    return {
        **source,
        "round_trip_digest": source_digest.hexdigest(),
        "sha256": file_sha256(path),
    }


def historical_classes(operator_name, historical_path):
    records = [
        item.force(trusted=True)
        for item in read_completions(historical_path, trusted=True)[operator_name]
    ]
    if operator_name in PROJECTED_LORENTZ_OPERATORS:
        basis = HistoricalDerivativeBasis.from_operator(
            DERIV_EFF_OPERATORS[operator_name]
        )
        for completion in records:
            completion.lorentz_projection = basis.project_existing_local(
                completion.operator.operator
            )
    classes = []
    append_unique_completions(classes, records)
    return records, classes


def _write_models(path, models):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
        delete=False,
    )
    temporary_path = Path(temporary.name)
    try:
        with temporary:
            for fields, topologies in sorted(models.items()):
                temporary.write(
                    json.dumps(
                        {
                            "fields": list(fields),
                            "topologies": [list(item) for item in sorted(topologies)],
                        },
                        sort_keys=True,
                        separators=(",", ":"),
                    )
                    + "\n"
                )
        temporary_path.replace(path)
    finally:
        temporary_path.unlink(missing_ok=True)
    return file_sha256(path)


def audit_exact_artifact(
    exact_path,
    model_path,
    historical,
    *,
    work_dir=None,
):
    """Audit exact classes and historical coverage with a disk-backed index."""

    exact_path = Path(exact_path)
    if work_dir is not None:
        Path(work_dir).mkdir(parents=True, exist_ok=True)
    stats = _new_stats()
    digest = OrderedSignatureDigest()
    models = defaultdict(set)
    comparisons = 0
    missing = []

    with tempfile.TemporaryDirectory(dir=work_dir) as directory:
        database_path = Path(directory) / "historical-coverage.sqlite3"
        with sqlite3.connect(database_path) as connection:
            connection.execute("PRAGMA temp_store = FILE")
            connection.execute("PRAGMA cache_size = -8192")
            connection.execute(
                "CREATE TABLE exact_classes (bucket_key TEXT, payload TEXT)"
            )
            for completion in iter_completion_jsonl(exact_path):
                validate_completion(completion)
                _update_stats(stats, completion)
                digest.update(round_trip_signature(completion))
                models[model_strings(completion)].add(
                    (completion.topology, completion.canonical_topology)
                )
                connection.execute(
                    "INSERT INTO exact_classes(bucket_key, payload) VALUES (?, ?)",
                    (
                        repr(exact_completion_bucket_key(completion)),
                        dumps_completion(completion),
                    ),
                )
            connection.execute(
                "CREATE INDEX exact_class_buckets ON exact_classes(bucket_key)"
            )
            connection.commit()

            for historical_completion in historical:
                clear_interaction_graph_cache()
                rows = connection.execute(
                    "SELECT payload FROM exact_classes WHERE bucket_key = ?",
                    (repr(exact_completion_bucket_key(historical_completion)),),
                )
                for (payload,) in rows:
                    comparisons += 1
                    if are_equivalent_completions(
                        historical_completion, loads_completion(payload)
                    ):
                        break
                else:
                    missing.append(completion_fingerprint(historical_completion))
                clear_interaction_graph_cache()

    model_sha256 = _write_models(model_path, models)
    return {
        **_serialise_stats(stats),
        "round_trip_digest": digest.hexdigest(),
        "democratic_models": len(models),
        "model_artifact": {
            "path": str(Path(model_path).resolve()),
            "sha256": model_sha256,
        },
        "historical_comparisons": comparisons,
        "missing_historical_fingerprints": [repr(item) for item in missing],
    }


def operator_scale_gev(operator):
    return float(
        max(
            numerical_np_scale_estimate(estimate)
            for estimate in neutrino_mass_estimate(operator)
        )
    )


def read_model_artifact(path):
    with Path(path).open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, start=1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
                fields = tuple(record["fields"])
                topologies = tuple(tuple(item) for item in record["topologies"])
            except (KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
                raise ValueError(
                    f"Invalid democratic-model record at {path}:{line_number}"
                ) from error
            if fields != tuple(sorted(set(fields))):
                raise ValueError(
                    f"Non-canonical democratic model at {path}:{line_number}"
                )
            yield fields, topologies


def generated_one_loop_weinberg_models():
    from neutrinomass.database.heavyloops import generate_models

    return tuple(
        sorted(
            {
                tuple(sorted({f"{field},0" for field in model}))
                for model in generate_models()
            }
        )
    )


def filter_democratic_registry(registry, scales, one_loop_models):
    """Filter a regenerated democratic registry by mass and one-loop Weinberg."""

    if set(registry) != set(scales):
        raise ValueError("Model registry and operator-scale keys differ")
    survivors = {
        operator: dict(models) for operator, models in registry.items()
    }
    removed_by_mass = Counter()

    def contains_registered_subset(fields, registered):
        fields = tuple(sorted(fields))
        return any(
            subset in registered
            for size in range(1, len(fields) + 1)
            for subset in combinations(fields, size)
        )

    active_sieves = set()
    scale_groups = defaultdict(list)
    for operator, scale in scales.items():
        scale_groups[scale].append(operator)
    for scale in sorted(scale_groups, reverse=True):
        group = sorted(scale_groups[scale])
        for operator in group:
            retained = {}
            for fields, topologies in survivors[operator].items():
                if contains_registered_subset(fields, active_sieves):
                    removed_by_mass[operator] += 1
                else:
                    retained[fields] = topologies
            survivors[operator] = retained
        for operator in group:
            active_sieves.update(survivors[operator])

    one_loop_scale = 605520000000.0 / (16 * math.pi**2)
    one_loop_sieves = set(one_loop_models)
    removed_by_one_loop = Counter()
    for operator in sorted(scales):
        if scales[operator] >= one_loop_scale:
            continue
        retained = {}
        for fields, topologies in survivors[operator].items():
            if contains_registered_subset(fields, one_loop_sieves):
                removed_by_one_loop[operator] += 1
            else:
                retained[fields] = topologies
        survivors[operator] = retained

    return {
        "survivors": survivors,
        "removed_by_mass": dict(sorted(removed_by_mass.items())),
        "removed_by_one_loop_weinberg": dict(
            sorted(removed_by_one_loop.items())
        ),
        "one_loop_scale_gev": one_loop_scale,
    }


def write_filtered_registry(path, survivors):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
        delete=False,
    )
    temporary_path = Path(temporary.name)
    count = 0
    try:
        with temporary:
            for operator, models in sorted(survivors.items()):
                for fields, topologies in sorted(models.items()):
                    count += 1
                    temporary.write(
                        json.dumps(
                            {
                                "operator": operator,
                                "fields": list(fields),
                                "topologies": [
                                    list(item) for item in sorted(topologies)
                                ],
                            },
                            sort_keys=True,
                            separators=(",", ":"),
                        )
                        + "\n"
                    )
        temporary_path.replace(path)
    finally:
        temporary_path.unlink(missing_ok=True)
    return {"records": count, "sha256": file_sha256(path)}


def census_operator(operator_name, historical_path, output_dir, *, hash_seed=None):
    """Generate and verify one registry operator from first principles."""

    started = perf_counter()
    registry = operator_registry()
    operator = registry[operator_name]
    kind = operator_kind(operator_name)
    derivative_count = (
        operator_strip_derivs(operator.operator)["n_derivs"]
        if kind == DERIVATIVE_KIND
        else 0
    )
    output_dir = Path(output_dir)
    historical_path = Path(historical_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_path = output_dir / f"op_{operator_name}_remediated.jsonl"
    exact_path = output_dir / f"op_{operator_name}_remediated_unique.jsonl"
    model_path = output_dir / f"op_{operator_name}_models.jsonl"

    generated = write_generated_artifact(
        raw_path, completion_stream(operator_name)
    )
    clear_interaction_graph_cache()
    deduplication = deduplicate_completion_jsonl(
        raw_path, exact_path, work_dir=output_dir
    )
    records, classes = historical_classes(operator_name, historical_path)
    exact = audit_exact_artifact(
        exact_path,
        model_path,
        classes,
        work_dir=output_dir,
    )
    if exact["missing_historical_fingerprints"]:
        raise ValueError(
            f"{len(exact['missing_historical_fingerprints'])} of "
            f"{len(classes)} historical classes missing"
        )
    if generated["records"] != deduplication["input_records"]:
        raise ValueError("deduplication input count differs from generated count")
    if exact["records"] != deduplication["exact_classes"]:
        raise ValueError("deduplication output count differs from exact audit")
    if generated["sha256"] != deduplication["source_sha256"]:
        raise ValueError("raw artifact hash changed before deduplication")
    exact_sha256 = file_sha256(exact_path)
    if exact_sha256 != deduplication["destination_sha256"]:
        raise ValueError("exact artifact hash changed after deduplication")

    return {
        "schema_version": 1,
        "operator": operator_name,
        "kind": kind,
        "derivatives": derivative_count,
        "hash_seed": str(
            hash_seed
            if hash_seed is not None
            else os.environ.get("PYTHONHASHSEED", "")
        ),
        "wall_time_seconds": perf_counter() - started,
        "peak_memory_mib": peak_memory_mib(),
        "operator_scale_gev": operator_scale_gev(operator),
        "records": {
            "generator": generated["records"],
            "local": generated["local"],
            "routed": generated["routed"],
        },
        "exact_classes": {
            "all": exact["records"],
            "local": exact["local"],
            "routed": exact["routed"],
        },
        "democratic_models": exact["democratic_models"],
        "historical": {
            "records": len(records),
            "classes": len(classes),
            "reproduced": len(classes),
            "missing": 0,
            "exact_comparisons": exact["historical_comparisons"],
            "artifact": {
                "path": str(historical_path.resolve()),
                "sha256": file_sha256(historical_path),
            },
        },
        "topologies": {
            "generator": generated["topologies"],
            "exact_classes": exact["topologies"],
        },
        "completion_digests": {
            "generator": deduplication["input_completion_digest"],
            "exact_classes": deduplication["exact_completion_digest"],
        },
        "round_trip_digests": {
            "generator": generated["round_trip_digest"],
            "exact_classes": exact["round_trip_digest"],
        },
        "deduplication": deduplication,
        "artifacts": {
            "generator": {
                "path": str(raw_path.resolve()),
                "sha256": generated["sha256"],
            },
            "exact_classes": {
                "path": str(exact_path.resolve()),
                "sha256": exact_sha256,
            },
            "democratic_models": exact["model_artifact"],
        },
    }
