#!/usr/bin/env python3
"""Reproduce the released model/filter contract without executing legacy tails."""

import argparse
import ast
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path

from neutrinomass.database.database import DATA, EXOTICS
from neutrinomass.database.filtering import filter_published_registry, published_model_key
from neutrinomass.database.rebuild import (
    operator_registry, operator_scale_tev, filter_democratic_registry,
    generated_one_loop_weinberg_models,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("legacy_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    registry, manifest, legacy_terms, literal_union = {}, {}, set(), set()
    first_terms={}
    for name in sorted(operator_registry()):
        path = args.legacy_dir / f"op_{name}.dat"
        models, records = {}, 0
        first_terms[name]={}
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for line in stream:
                digest.update(line)
                if not line.strip():
                    continue
                text = line.decode("utf-8")
                prefix = "LazyCompletion(head="
                if not text.startswith(prefix):
                    raise ValueError(f"Unexpected legacy record in {path}")
                head = ast.literal_eval(text[len(prefix):].split(", tail=", 1)[0])
                if head["operator_name"] != name:
                    raise ValueError("Legacy operator mismatch")
                fields = tuple(sorted(set(head["quantum_numbers"])))
                models[fields] = ()
                # The published stringent key forgets contractions and labels,
                # but retains repeated interaction factors and term conjugacy.
                from neutrinomass.utils.functions import conjugate_term
                terms = tuple(sorted(min(tuple(sorted(t)), conjugate_term(tuple(t)))
                                     for t in head["terms"]))
                legacy_terms.add(terms)
                first_terms[name].setdefault(fields,terms)
                records += 1
        registry[name] = models
        literal_union.update(models)
        manifest[name] = dict(sha256=digest.hexdigest(), records=records, models=len(models))
        print(json.dumps(dict(event="baseline_input", operator=name, records=records)), flush=True)
    operators = operator_registry()
    scales = {name: operator_scale_tev(op) for name, op in operators.items()}
    result = filter_published_registry(registry, scales)
    baseline = defaultdict(set)
    for operator, number in zip(DATA["op"], DATA["democratic_num"]):
        baseline[operator].add(int(number))
    mismatches = {}
    for operator, models in result["survivors"].items():
        actual = {math.prod(int(EXOTICS[field]) for field in fields) for fields in models}
        if actual != baseline[operator]:
            mismatches[operator] = dict(added=sorted(actual-baseline[operator]),
                                        removed=sorted(baseline[operator]-actual))
    legacy_v10 = filter_democratic_registry(registry, scales, generated_one_loop_weinberg_models())
    v10_union = set().union(*(set(m) for m in legacy_v10["survivors"].values()))
    comparison = {
        name: dict(
            v10_only=[list(f) for f in sorted(set(legacy_v10["survivors"][name])-set(result["survivors"][name]))],
            published_only=[list(f) for f in sorted(set(result["survivors"][name])-set(legacy_v10["survivors"][name]))])
        for name in registry
        if set(legacy_v10["survivors"][name]) != set(result["survivors"][name])
    }
    filtered_first_terms={first_terms[name][fields] for name,models in result["survivors"].items() for fields in models}
    report = dict(schema_version=2, scale_unit="TeV", inputs=manifest,
        literal_unfiltered_models=len(literal_union),
        conjugacy_unfiltered_models=result["distinct_input_models"],
        filtered_models=result["distinct_surviving_models"],
        per_operator_differences=mismatches, v10_filtered_models=len(v10_union),
        v10_filter_differences=comparison, historical_interaction_multisets=len(legacy_terms),
        published_lagrangians=430810,
        historical_filtered_first_representative_lagrangians=len(filtered_first_terms),
        packaged_filtered_lagrangians=len(set(DATA["stringent_num"])),
        historical_interaction_scope="243 raw archive operators; contraction-blind published term keys",
        explanations=["Published conjugate fields share identity.",
            "Published loop prefix lookup inherits archived baryon assignments; v10 forces 3B=0.",
            "Published loop templates retain repeated factors; v10 collapses them."])
    report["passed"] = (report["literal_unfiltered_models"] == 141989
        and report["filtered_models"] == 11216 and not mismatches
        and report["historical_interaction_multisets"] == 430810
        and report["historical_filtered_first_representative_lagrangians"] == 11483
        and not (set(baseline)-set(registry)))
    (args.output_dir/"baseline.json").write_text(json.dumps(report, indent=2, sort_keys=True)+"\n")
    with (args.output_dir/"published-models.jsonl").open("w") as stream:
        for name, models in sorted(result["survivors"].items()):
            for fields in sorted(models):
                stream.write(json.dumps(dict(operator=name, fields=fields,
                    key=published_model_key(fields)), sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in report.items() if k not in {"inputs","v10_filter_differences"}}))
    if not report["passed"]:
        raise SystemExit("Published baseline did not reproduce")


if __name__ == "__main__":
    main()
