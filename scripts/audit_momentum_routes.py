"""Read-only endpoint audit. Run with the intended checkout on PYTHONPATH."""
import argparse
from collections import Counter
import hashlib
import gzip
import json
from pathlib import Path
import subprocess
import sys

import neutrinomass
from neutrinomass.database.serialization import loads_completion
from neutrinomass.completions.routing import inspect_endpoints, migrate_legacy_one_derivative
from neutrinomass.completions.completions import canonical_propagator_cut, partition_leaves


def current(completion, edge):
    leaves = {leaf.node: leaf.field for leaf in partition_leaves(completion.partition)}
    value = 0
    for node in canonical_propagator_cut(completion.partition, completion.graph, edge):
        field = leaves[node]
        field = field.strip_derivs() if field.derivs else field
        if field.label.rstrip("†") == "phidm":
            value += -1 if field.is_conj else 1
    return value != 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("streams", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = {"version": 1, "validator": "physical-endpoints-v1-development",
              "interpreter": sys.executable, "package": neutrinomass.__file__,
              "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
              "dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], text=True)),
              "streams": []}
    for path in args.streams:
        raw = path.read_bytes()
        counts = Counter()
        records = []
        decoded = gzip.decompress(raw) if path.suffix == ".gz" else raw
        for line, payload in enumerate(decoded.splitlines(), 1):
            counts["raw"] += 1
            try:
                completion = loads_completion(payload)
                audit = inspect_endpoints(completion)
                counts["consistent" if audit.agrees else "inconsistent"] += 1
                row = {"line": line, "physical": audit.momentum_edges,
                       "stored": audit.stored_momentum_edges,
                       "status": "consistent" if audit.agrees else "invalid_stored_route"}
                if len(audit.momentum_edges) == len(audit.stored_momentum_edges) == 1:
                    before = current(completion, audit.stored_momentum_edges[0])
                    after = current(completion, audit.momentum_edges[0])
                    counts[f"current_{int(before)}_to_{int(after)}"] += 1
                    migrate_legacy_one_derivative(completion)
                    counts["raw_migratable"] += 1
                    row.update(current_before=before, current_after=after,
                               raw_migration="validated_graph_only")
                elif not audit.agrees:
                    row["status"] = "unsupported_missing_or_multiple_routes"
                records.append(row)
            except Exception as error:
                counts["unsupported"] += 1
                records.append({"line": line, "unsupported": f"{type(error).__name__}: {error}"})
        manifest = path.with_suffix(".manifest.json")
        report["streams"].append({"path": str(path.resolve()),
            "sha256": hashlib.sha256(raw).hexdigest(), "counts": dict(counts),
            "producer_manifest": json.loads(manifest.read_text()) if manifest.exists() else None,
            "records": records})
        print(path.name, dict(counts), flush=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
