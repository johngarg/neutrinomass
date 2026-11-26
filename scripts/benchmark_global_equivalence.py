"""Compare exact global matching with a frozen Git revision on saved pairs.

Pairs are JSON arrays of two serialized completion strings. Every negative
pair also supplies a positive control with relabelled species and reversed
interaction order. Timings exclude decoding and use independent graph caches.
The reference implementation is saved alongside the result for replay.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import time

from neutrinomass.completions import equivalence
from neutrinomass.database import global_equivalence
from neutrinomass.database.checkpoints import atomic_json, digest_file
from neutrinomass.database.provenance import runtime_provenance
from neutrinomass.database.serialization import dumps_completion, loads_completion


def reference_module(root, revision, relative, destination):
    source = subprocess.check_output(["git", "show", f"{revision}:{relative}"], cwd=root)
    destination.write_bytes(source)
    spec = importlib.util.spec_from_file_location(destination.stem, destination)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def relabelled(completion):
    labels = sorted(global_equivalence.species_descriptors(completion))
    mapping = {label: f"benchmark_species_{i}" for i, label in enumerate(labels)}

    def transform(value):
        if isinstance(value, dict):
            result = {k: transform(v) for k, v in value.items()}
            if isinstance(result.get("label"), str):
                result["label"] = equivalence.remapped_field_label(result["label"], mapping)
            return result
        if isinstance(value, list):
            return [transform(v) for v in value]
        return value

    result = loads_completion(json.dumps(transform(json.loads(dumps_completion(completion)))))
    result.terms = list(reversed(result.terms))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pairs", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--reference-revision", required=True)
    args = parser.parse_args()
    provenance = runtime_provenance()
    root = Path(provenance["source_directory"])
    revision = subprocess.check_output(
        ["git", "rev-parse", args.reference_revision + "^{commit}"], cwd=root, text=True).strip()
    args.output.mkdir(parents=True, exist_ok=False)
    ref_graph = reference_module(root, revision, "neutrinomass/completions/equivalence.py",
                                 args.output / "reference_equivalence.py")
    ref_global = reference_module(root, revision, "neutrinomass/database/global_equivalence.py",
                                  args.output / "reference_global_equivalence.py")
    ref_global.equivalent_lagrangians = ref_graph.equivalent_lagrangians
    pairs = [tuple(map(loads_completion, pair)) for pair in json.loads(args.pairs.read_text())]
    cases = pairs + [(left, relabelled(left)) for left, _ in pairs]
    runs = {}
    for name, module, graph_module in (
            ("reference", ref_global, ref_graph),
            ("optimised", global_equivalence, equivalence)):
        graph_module.clear_interaction_graph_cache()
        start, cpu = time.monotonic(), time.process_time()
        answers = [module.equivalent_global_lagrangians(a, b) for a, b in cases]
        runs[name] = dict(wall_seconds=time.monotonic()-start,
                          cpu_seconds=time.process_time()-cpu, answers=answers)
    agreement = runs["reference"]["answers"] == runs["optimised"]["answers"]
    positives = all(runs["optimised"]["answers"][len(pairs):])
    summary = dict(complete=agreement and positives, agreement=agreement,
                   positive_controls_passed=positives, cases=len(cases),
                   pairs_sha256=digest_file(args.pairs), reference_revision=revision,
                   producer=provenance, runs=runs)
    atomic_json(args.output / "benchmark.json", summary)
    print(json.dumps(summary, sort_keys=True), flush=True)
    if not summary["complete"]:
        raise AssertionError("Global-equivalence benchmark disagrees")


if __name__ == "__main__":
    main()
