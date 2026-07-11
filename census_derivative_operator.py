#!/usr/bin/env python3

"""Reproducible census for a supported one-derivative operator."""

import argparse
from collections import Counter, defaultdict
from functools import reduce
from hashlib import sha256
import json
from math import pi
from operator import mul
from pathlib import Path
import resource
import sys
from time import perf_counter

from neutrinomass.completions.completions import (
    append_unique_completions,
    are_equivalent_completions,
    base_exotic_label,
    deriv_operator_completions,
    is_singlet,
)
from neutrinomass.completions.fingerprints import (
    completion_digest,
    completion_fingerprint,
    democratic_model_fingerprint,
)
from neutrinomass.completions.operators import DERIV_EFF_OPERATORS
from neutrinomass.database import (
    MVDF,
    neutrino_mass_estimate,
    numerical_np_scale_estimate,
    read_completion_jsonl,
    read_completions,
    write_completion_jsonl,
)
from neutrinomass.database.heavyloops import generate_models


def file_digest(path):
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


def quantum_number_string(info):
    lorentz, colour_up, colour_down, isospin, (_, baryon), (_, hypercharge) = info
    return (
        f"{lorentz},{colour_up}{colour_down},{isospin},"
        f"{hypercharge},{baryon}"
    )


def model_strings(completion):
    return tuple(
        sorted(
            quantum_number_string(info)
            for info in democratic_model_fingerprint(completion)
        )
    )


def model_number(model):
    missing = [field for field in model if field not in MVDF.exotics]
    if missing:
        raise ValueError(f"fields absent from packaged prime registry: {missing}")
    return reduce(mul, (MVDF.exotics[field] for field in model), 1)


def generated_one_loop_weinberg_numbers():
    numbers = set()
    for model in generate_models():
        fields = tuple(f"{field},0" for field in model)
        if all(field in MVDF.exotics for field in fields):
            numbers.add(reduce(mul, {MVDF.exotics[field] for field in fields}, 1))
    return numbers


def topology_distribution(completions):
    counts = Counter(
        f"{completion.topology}|{completion.canonical_topology}"
        for completion in completions
    )
    return dict(sorted(counts.items()))


def validate_physics(completions):
    for completion in completions:
        if any(term.safe_simplify() == 0 for term in completion.terms):
            raise ValueError("vanishing UV interaction")
        if any(
            not is_singlet(term)
            or sum(field.mass_dim for field in term.fields) > 4
            for term in completion.terms
        ):
            raise ValueError("non-renormalisable or non-singlet UV interaction")

        if not completion.derivative_routes:
            continue
        if len(completion.derivative_routes) != 1:
            raise ValueError("routed completion does not have exactly one route")
        route = completion.derivative_routes[0]
        if not completion.graph.has_edge(*route.edge):
            raise ValueError(f"recorded route edge is absent: {route.edge}")
        particle = completion.graph.edges[route.edge]["particle"]
        if base_exotic_label(particle) != base_exotic_label(route.numerator_field):
            raise ValueError("route numerator does not match the edge particle")


def verify_round_trip(path, expected):
    restored = read_completion_jsonl(path)
    if [completion_fingerprint(item) for item in restored] != [
        completion_fingerprint(item) for item in expected
    ]:
        raise ValueError(f"fingerprint changed across JSONL round trip: {path}")
    if any(
        not are_equivalent_completions(left, right)
        for left, right in zip(expected, restored)
    ):
        raise ValueError(f"equivalence changed across JSONL round trip: {path}")
    if any(
        left.derivative_routes != right.derivative_routes
        for left, right in zip(expected, restored)
    ):
        raise ValueError(f"route metadata changed across JSONL round trip: {path}")
    return restored


def filter_models(operator, completions):
    scale = max(
        numerical_np_scale_estimate(estimate)
        for estimate in neutrino_mass_estimate(operator)
    )
    upstream = {
        int(number)
        for number in MVDF.loc[MVDF["scale"] > scale, "democratic_num"]
    }
    one_loop = generated_one_loop_weinberg_numbers()

    topologies = defaultdict(set)
    for completion in completions:
        topologies[model_strings(completion)].add(
            (completion.topology, completion.canonical_topology)
        )

    survivors = []
    for model, model_topologies in sorted(topologies.items()):
        number = model_number(model)
        if any(number % sieve == 0 for sieve in upstream):
            continue
        if scale < 605520000000.0 / (16 * pi**2) and any(
            number % sieve == 0 for sieve in one_loop
        ):
            continue
        survivors.append(
            {
                "fields": list(model),
                "topologies": [list(item) for item in sorted(model_topologies)],
            }
        )

    return float(scale), survivors


def run(operator_name, historical_path, output_dir):
    started = perf_counter()
    operator = DERIV_EFF_OPERATORS[operator_name]
    raw = deriv_operator_completions(operator)
    validate_physics(raw)

    unique = []
    append_unique_completions(unique, raw)

    historical_records = [
        item.force(trusted=True)
        for item in read_completions(historical_path, trusted=True)[operator_name]
    ]
    historical_classes = []
    append_unique_completions(historical_classes, historical_records)
    missing = [
        completion
        for completion in historical_classes
        if not any(
            are_equivalent_completions(completion, candidate)
            for candidate in unique
        )
    ]
    if missing:
        raise ValueError(
            f"{len(missing)} of {len(historical_classes)} historical classes missing"
        )

    output_dir = Path(output_dir)
    raw_path = output_dir / f"op_{operator_name}_remediated.jsonl"
    unique_path = output_dir / f"op_{operator_name}_remediated_unique.jsonl"
    write_completion_jsonl(raw_path, raw)
    write_completion_jsonl(unique_path, unique)
    restored_raw = verify_round_trip(raw_path, raw)
    restored_unique = verify_round_trip(unique_path, unique)

    scale, survivors = filter_models(operator, unique)
    routed_raw = [item for item in raw if item.derivative_routes]
    routed_unique = [item for item in unique if item.derivative_routes]
    report = {
        "operator": operator_name,
        "hash_seed": __import__("os").environ.get("PYTHONHASHSEED"),
        "wall_time_seconds": perf_counter() - started,
        "peak_memory_mib": peak_memory_mib(),
        "records": {
            "generator": len(raw),
            "local": len(raw) - len(routed_raw),
            "routed": len(routed_raw),
        },
        "exact_classes": {
            "all": len(unique),
            "local": len(unique) - len(routed_unique),
            "routed": len(routed_unique),
        },
        "democratic_models": len({model_strings(item) for item in unique}),
        "historical": {
            "records": len(historical_records),
            "classes": len(historical_classes),
            "reproduced": len(historical_classes),
            "missing": 0,
        },
        "topologies": {
            "generator": topology_distribution(raw),
            "exact_classes": topology_distribution(unique),
        },
        "filtering": {
            "operator_scale_gev": scale,
            "survivor_count": len(survivors),
            "survivors": survivors,
        },
        "completion_digests": {
            "generator": completion_digest(restored_raw),
            "exact_classes": completion_digest(restored_unique),
        },
        "artifacts": {
            "generator": {
                "path": str(raw_path.resolve()),
                "sha256": file_digest(raw_path),
            },
            "exact_classes": {
                "path": str(unique_path.resolve()),
                "sha256": file_digest(unique_path),
            },
        },
    }
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operator", choices=sorted(DERIV_EFF_OPERATORS))
    parser.add_argument("historical_path", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()

    report = run(args.operator, args.historical_path, args.output_dir)
    payload = json.dumps(report, indent=2, sort_keys=True)
    if args.report:
        args.report.write_text(payload + "\n", encoding="utf-8")
    print(payload)


if __name__ == "__main__":
    main()
