"""Explain changed record counts by exact migrated-graph comparisons."""
from collections import Counter
from copy import copy
import json
from pathlib import Path
import sys

from neutrinomass.database.serialization import completion_from_record
from neutrinomass.completions.routing import inspect_endpoints, migrate_legacy_one_derivative
from neutrinomass.completions.completions import are_equivalent_completions, append_unique_completions
from neutrinomass.completions.fingerprints import completion_digest

directory = Path(sys.argv[1])
before = json.loads((directory / "benchmarks-before.json").read_text())
after = json.loads((directory / "benchmarks-after.json").read_text())
result = {}
for operator in ["D6a", "D3"]:
    old = [completion_from_record(r) for r in before[operator]]
    new = [completion_from_record(r) for r in after[operator]]
    rows = []
    for index, completion in enumerate(old):
        if not completion.derivative_routes:
            continue
        audit = inspect_endpoints(completion)
        migrated = migrate_legacy_one_derivative(completion)
        matches = []
        for target, candidate in enumerate(new):
            candidate = copy(candidate)
            candidate.lorentz_projection = None
            if are_equivalent_completions(migrated, candidate):
                matches.append(target)
        rows.append({"old_index": index, "agrees": audit.agrees,
            "stored": audit.stored_momentum_edges, "physical": audit.momentum_edges,
            "migrated_graph_matches": matches})
    unique = []
    append_unique_completions(unique, new)
    result[operator] = {"old_count": len(old), "new_count": len(new),
        "routed": sum(bool(c.derivative_routes) for c in new),
        "unique_count": len(unique), "unique_routed": sum(bool(c.derivative_routes) for c in unique),
        "topologies": dict(Counter(c.topology for c in new)),
        "digest": completion_digest(new), "unique_digest": completion_digest(unique),
        "old_route_audit": rows}
(directory / "benchmark-comparison.json").write_text(json.dumps(result, indent=2) + "\n")
