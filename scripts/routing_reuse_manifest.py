"""Conservative local artifact inventory for the erratum census handoff."""
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys

import neutrinomass

checkout = Path(__file__).resolve().parents[1]
vault = checkout.parents[2]
output = vault / ".Data/exploding-ops/routing-20261030/reuse-manifest.json"
files = list(checkout.parent.glob("op_D*"))
files += list((vault / ".Data/exploding-ops/v10-review-20260913/imported-v10").glob("seed-*/operators/D*/*"))
inventory = []
for path in sorted(set(files)):
    if not path.is_file():
        continue
    inventory.append({"path": str(path.relative_to(vault)),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "action": "retain_as_historical_evidence_not_production_reuse",
        "reason": "Pre-repair derivative output; upstream acceptance/equivalence may have discarded candidates. Full raw-stage coverage not certified."})
audits = []
for name in ["dm-endpoints.json", "legacy-benchmarks.json"]:
    for stream in json.loads((output.parent / name).read_text())["streams"]:
        counts = stream["counts"]
        audits.append({"path": stream["path"], "old_sha256": stream["sha256"],
            "producer_manifest": stream["producer_manifest"], "counts": counts,
            "validator": "physical-endpoints-v1-development",
            "action": "raw_graph_migration_only" if counts.get("raw_migratable") == counts["raw"]
                      else "requires_generation_or_explicit_unsupported_classification",
            "downstream": ["projection", "route-sensitive equivalence", "amplitude audit", "filtering", "published comparison"],
            "coverage": "not_certified_by_record_repair"})
sources = {}
for path in sorted(checkout.glob("neutrinomass/**/*.py")):
    sources[str(path.relative_to(checkout))] = hashlib.sha256(path.read_bytes()).hexdigest()
manifest = {"schema": "neutrinomass.routing-reuse", "version": 1,
    "production_reuse_approved": False,
    "base_commit": "d77831f645d926b3502bd71c76ad04050d7c20dd",
    "validator_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=checkout, text=True).strip(),
    "branch": subprocess.check_output(["git", "branch", "--show-current"], cwd=checkout, text=True).strip(),
    "interpreter": sys.executable, "python": platform.python_version(), "package": neutrinomass.__file__,
    "requirements_lock_sha256": hashlib.sha256((checkout / "requirements-test.lock").read_bytes()).hexdigest(),
    "source_sha256": sources, "audited_streams": audits, "local_derivative_artifacts": inventory,
    "schema_policy": {"write": 3, "read": [1, 2, 3],
        "legacy_read": "Decode for audit only; no implicit route repair or retained projection certification.",
        "migration": "migrate_legacy_one_derivative requires complete unambiguous endpoints and a single supplied momentum contribution; clears Lorentz projection; retain original bytes and producer.",
        "limits": "Version 3 validates supplied routes, not completeness of legacy local-derivative enumeration. No global certification."},
    "remaining": ["Local derivative completions with physical momentum but absent route metadata",
        "Full D20 and D5b regenerated benchmarks and historical coverage",
        "Complete supported multi-derivative generation validation",
        "Full derivative operator/seed affected-stage inventory including remote artifacts",
        "DM catalogue regeneration and combined weak-projector/class comparisons"],
    "scope": "Local conservative handoff; remote registry not inspected; no production jobs submitted."}
output.write_text(json.dumps(manifest, indent=2) + "\n")
print(output, len(inventory), "local artifacts")
