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
    UnreducedSecondDerivativeBasis,
    append_unique_completions,
    are_equivalent_completions,
    base_exotic_label,
    canonical_propagator_cut,
    deriv_operator_completions,
    exact_completion_bucket_key,
    exotic_species,
    is_singlet,
    operator_completions,
    operator_strip_derivs,
    unique_multi_derivative_projection,
)
from neutrinomass.completions.amplitudes import (
    audit_amplitude_symmetrisation,
)
from neutrinomass.completions.equivalence import clear_interaction_graph_cache
from neutrinomass.completions.fingerprints import (
    completion_fingerprint,
    democratic_model_fingerprint,
    propagator_model_fingerprint,
    species_model_fingerprint,
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
        operator = DERIV_EFF_OPERATORS[operator_name]
        derivative_count = operator_strip_derivs(operator.operator)["n_derivs"]
        return iter(
            deriv_operator_completions(
                operator,
                canonical_partitions=derivative_count <= 1,
            )
        )
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


def validate_completion(completion, *, check_vanishing=True):
    exotic_species(completion)
    if check_vanishing and any(
        term.safe_simplify() == 0 for term in completion.terms
    ):
        raise ValueError("vanishing UV interaction")
    if any(
        not is_singlet(term)
        or sum(field.mass_dim for field in term.fields) > 4
        for term in completion.terms
    ):
        raise ValueError("non-renormalisable or non-singlet UV interaction")

    contributions = getattr(completion, "momentum_contributions", ())
    if not contributions:
        return
    for contribution in contributions:
        if contribution.derivative_degree <= 0:
            raise ValueError("recorded propagator contribution has degree zero")
        if not completion.graph.has_edge(*contribution.edge):
            raise ValueError(
                f"recorded contribution edge is absent: {contribution.edge}"
            )
        particle = completion.graph.edges[contribution.edge]["particle"]
        if base_exotic_label(particle) != base_exotic_label(
            contribution.particle
        ):
            raise ValueError(
                "propagator contribution does not match the edge particle"
            )
        exact_matches = [
            field
            for field in completion.exotics
            if field.label == contribution.particle
        ]
        matching_fields = exact_matches or [
            field
            for field in completion.exotics
            if base_exotic_label(field.label)
            == base_exotic_label(contribution.particle)
        ]
        if (
            not matching_fields
            or len({field.is_fermion for field in matching_fields}) != 1
        ):
            raise ValueError("propagator contribution has ambiguous particle")
        field = matching_fields[0]
        if contribution.numerator_kind == "scalar" and not field.is_boson:
            raise ValueError("scalar numerator recorded on a fermion edge")
        if (
            contribution.numerator_kind in {"mass", "momentum"}
            and not field.is_fermion
        ):
            raise ValueError("fermion numerator recorded on a scalar edge")
        if (
            contribution.cut_side
            and contribution.cut_side
            != canonical_propagator_cut(
                completion.partition, completion.graph, contribution.edge
            )
        ):
            raise ValueError("propagator contribution has a noncanonical cut")

    operator = DERIV_EFF_OPERATORS.get(completion.operator.name)
    if operator is not None:
        requested_degree = operator_strip_derivs(operator.operator)["n_derivs"]
        generated_degree = sum(
            contribution.derivative_degree for contribution in contributions
        )
        if generated_degree != requested_degree:
            raise ValueError(
                "propagator contribution degree does not match the operator"
            )


def round_trip_signature(completion):
    projection = getattr(completion, "lorentz_projection", None)
    return (
        completion_fingerprint(completion),
        tuple(completion.momentum_contributions),
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
    routed = bool(completion.momentum_contributions)
    stats["routed" if routed else "local"] += 1
    stats["topologies"][topology_key(completion)] += 1


def _serialise_stats(stats):
    return {
        "records": stats["records"],
        "local": stats["local"],
        "routed": stats["routed"],
        "topologies": dict(sorted(stats["topologies"].items())),
    }


def write_generated_artifact(path, completions, *, terms_prevalidated=False):
    """Atomically write, validate and stream-audit generated completions.

    ``terms_prevalidated`` is reserved for ``completion_stream``: the partition
    constructor has already rejected every symbolically vanishing vertex.
    Singlet, mass-dimension and route invariants are always checked here.
    """

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
                validate_completion(
                    completion, check_vanishing=not terms_prevalidated
                )
                _update_stats(source_stats, completion)
                source_digest.update(round_trip_signature(completion))
                temporary.write(dumps_completion(completion) + "\n")
        temporary_path.replace(path)
    finally:
        temporary_path.unlink(missing_ok=True)

    decoded_stats = _new_stats()
    decoded_digest = OrderedSignatureDigest()
    # Source objects were fully validated before writing.  This second pass
    # isolates schema/metadata round-trip checks without repeating symbolic
    # simplification for every UV term.
    for completion in iter_completion_jsonl(path):
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


def _invalid_historical_class(completion, error):
    return {
        "fingerprint": repr(completion_fingerprint(completion)),
        "topology": topology_key(completion),
        "reason": str(error),
        "vanishing_term_indices": [
            index
            for index, term in enumerate(completion.terms)
            if term.safe_simplify() == 0
        ],
    }


def partition_bucketable_historical_records(records):
    """Separate records that cannot enter physical-equivalence bucketing."""

    bucketable = []
    invalid = []
    for completion in records:
        try:
            exact_completion_bucket_key(completion)
        except ValueError:
            try:
                validate_completion(completion)
            except ValueError as error:
                invalid.append(_invalid_historical_class(completion, error))
            else:
                raise
        else:
            bucketable.append(completion)
    return bucketable, invalid


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
    else:
        derivative_operator = DERIV_EFF_OPERATORS.get(operator_name)
        projection = (
            unique_multi_derivative_projection(derivative_operator)
            if derivative_operator is not None
            else None
        )
        if projection is not None:
            for completion in records:
                completion.lorentz_projection = projection
        elif derivative_operator is not None and operator_strip_derivs(
            derivative_operator.operator
        )["n_derivs"] == 2:
            basis = UnreducedSecondDerivativeBasis.from_operator(
                derivative_operator
            )
            for completion in records:
                completion.lorentz_projection = basis.project_existing_local(
                    completion.operator.operator
                )
    bucketable, invalid = partition_bucketable_historical_records(records)
    classes = []
    append_unique_completions(classes, bucketable)
    return records, classes, invalid


def classify_historical_classes(classes):
    valid = []
    invalid = []
    for completion in classes:
        try:
            validate_completion(completion)
        except ValueError as error:
            invalid.append(_invalid_historical_class(completion, error))
        else:
            audit = audit_amplitude_symmetrisation(completion)
            if audit.status == "unsupported":
                raise ValueError(
                    "Cannot audit historical completion amplitude: "
                    f"{audit.reason}"
                )
            if audit.is_zero:
                invalid.append(
                    _invalid_historical_class(completion, ValueError(audit.reason))
                )
            else:
                valid.append(completion)
    return valid, invalid


def audit_amplitude_artifact(source, destination):
    """Filter structurally exact classes by full amplitude symmetrisation."""

    source = Path(source)
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        prefix=f".{destination.name}.",
        suffix=".tmp",
        dir=destination.parent,
        delete=False,
    )
    temporary_path = Path(temporary.name)
    fingerprint_database = tempfile.NamedTemporaryFile(
        prefix=f".{destination.name}.fingerprints.",
        suffix=".sqlite3",
        dir=destination.parent,
        delete=False,
    )
    fingerprint_database_path = Path(fingerprint_database.name)
    fingerprint_database.close()
    input_records = 0
    output_records = 0
    rejected_topologies = Counter()
    derivative_sectors = Counter()
    try:
        with temporary, sqlite3.connect(fingerprint_database_path) as connection:
            connection.execute("CREATE TABLE fingerprints (value TEXT NOT NULL)")
            for completion in iter_completion_jsonl(source):
                input_records += 1
                audit = audit_amplitude_symmetrisation(completion)
                if audit.status == "unsupported":
                    raise ValueError(
                        "Cannot audit completion amplitude at record "
                        f"{input_records}: {audit.reason}"
                    )
                derivative_sectors[audit.derivative_sectors] += 1
                if audit.is_zero:
                    rejected_topologies[topology_key(completion)] += 1
                    continue
                output_records += 1
                temporary.write(dumps_completion(completion) + "\n")
                connection.execute(
                    "INSERT INTO fingerprints(value) VALUES (?)",
                    (repr(completion_fingerprint(completion)),),
                )
            connection.commit()
            completion_digest = sha256()
            first = True
            for (fingerprint,) in connection.execute(
                "SELECT value FROM fingerprints ORDER BY value"
            ):
                if not first:
                    completion_digest.update(b"\n")
                completion_digest.update(fingerprint.encode("utf-8"))
                first = False
        temporary_path.replace(destination)
    finally:
        temporary_path.unlink(missing_ok=True)
        fingerprint_database_path.unlink(missing_ok=True)

    return {
        "input_records": input_records,
        "surviving_records": output_records,
        "rejected_records": input_records - output_records,
        "rejected_topologies": dict(sorted(rejected_topologies.items())),
        "derivative_sectors": {
            str(key): value for key, value in sorted(derivative_sectors.items())
        },
        "completion_digest": completion_digest.hexdigest(),
        "source_sha256": file_sha256(source),
        "destination_sha256": file_sha256(destination),
    }


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
    species_models = set()
    propagator_models = set()
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
            # Deduplication only removes records from the fully validated raw
            # stream, so the exact audit need only verify structure, metadata,
            # models and historical coverage.
            for completion in iter_completion_jsonl(exact_path):
                _update_stats(stats, completion)
                digest.update(round_trip_signature(completion))
                models[model_strings(completion)].add(
                    (completion.topology, completion.canonical_topology)
                )
                species_models.add(species_model_fingerprint(completion))
                propagator_models.add(
                    propagator_model_fingerprint(completion)
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
        "species_models": len(species_models),
        "propagator_models": len(propagator_models),
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
    structural_path = (
        output_dir / f"op_{operator_name}_remediated_structural_unique.jsonl"
    )
    exact_path = output_dir / f"op_{operator_name}_remediated_unique.jsonl"
    model_path = output_dir / f"op_{operator_name}_models.jsonl"

    generated = write_generated_artifact(
        raw_path,
        completion_stream(operator_name),
        terms_prevalidated=True,
    )
    clear_interaction_graph_cache()
    deduplication = deduplicate_completion_jsonl(
        raw_path, structural_path, work_dir=output_dir
    )
    amplitude_audit = audit_amplitude_artifact(structural_path, exact_path)
    records, classes, unbucketable_historical = historical_classes(
        operator_name, historical_path
    )
    valid_historical, invalid_historical = classify_historical_classes(classes)
    invalid_historical = unbucketable_historical + invalid_historical
    exact = audit_exact_artifact(
        exact_path,
        model_path,
        valid_historical,
        work_dir=output_dir,
    )
    if exact["missing_historical_fingerprints"]:
        raise ValueError(
            f"{len(exact['missing_historical_fingerprints'])} of "
            f"{len(valid_historical)} valid historical classes missing"
        )
    if generated["records"] != deduplication["input_records"]:
        raise ValueError("deduplication input count differs from generated count")
    if amplitude_audit["input_records"] != deduplication["exact_classes"]:
        raise ValueError("amplitude-audit input differs from structural exact count")
    if exact["records"] != amplitude_audit["surviving_records"]:
        raise ValueError("amplitude-audit output count differs from exact audit")
    if generated["sha256"] != deduplication["source_sha256"]:
        raise ValueError("raw artifact hash changed before deduplication")
    structural_sha256 = file_sha256(structural_path)
    if structural_sha256 != deduplication["destination_sha256"]:
        raise ValueError("structural artifact hash changed after deduplication")
    exact_sha256 = file_sha256(exact_path)
    if exact_sha256 != amplitude_audit["destination_sha256"]:
        raise ValueError("physical exact artifact hash changed after amplitude audit")

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
        "structural_exact_classes": deduplication["exact_classes"],
        "democratic_models": exact["democratic_models"],
        "species_models": exact["species_models"],
        "propagator_models": exact["propagator_models"],
        "historical": {
            "records": len(records),
            "classes": len(classes) + len(unbucketable_historical),
            "valid_classes": len(valid_historical),
            "invalid_classes": len(invalid_historical),
            "invalid": invalid_historical,
            "reproduced": len(valid_historical),
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
            "structural_exact_classes": deduplication[
                "exact_completion_digest"
            ],
            "exact_classes": amplitude_audit["completion_digest"],
        },
        "round_trip_digests": {
            "generator": generated["round_trip_digest"],
            "exact_classes": exact["round_trip_digest"],
        },
        "deduplication": deduplication,
        "amplitude_symmetrisation": amplitude_audit,
        "artifacts": {
            "generator": {
                "path": str(raw_path.resolve()),
                "sha256": generated["sha256"],
            },
            "exact_classes": {
                "path": str(exact_path.resolve()),
                "sha256": exact_sha256,
            },
            "structural_exact_classes": {
                "path": str(structural_path.resolve()),
                "sha256": structural_sha256,
            },
            "democratic_models": exact["model_artifact"],
        },
    }
