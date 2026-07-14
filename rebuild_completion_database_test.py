import json

import pytest

from neutrinomass.database.rebuild import file_sha256
from rebuild_completion_database import (
    build_manifest,
    completed_report_is_valid,
    stable_report_data,
)


def sample_report(directory, seed):
    directory.mkdir(parents=True)
    artifacts = {}
    for name in (
        "generator",
        "structural_exact_classes",
        "exact_classes",
        "democratic_models",
    ):
        path = directory / f"{name}.jsonl"
        path.write_text(f"{name}\n", encoding="utf-8")
        artifacts[name] = {
            "path": str(path.resolve()),
            "sha256": file_sha256(path),
        }
    return {
        "schema_version": 1,
        "operator": "1",
        "kind": "regular",
        "derivatives": 0,
        "hash_seed": str(seed),
        "wall_time_seconds": 1.0 + seed,
        "peak_memory_mib": 2.0 + seed,
        "operator_scale_gev": 3.0,
        "records": {"generator": 8, "local": 8, "routed": 0},
        "structural_exact_classes": 4,
        "exact_classes": {"all": 3, "local": 3, "routed": 0},
        "amplitude_symmetrisation": {
            "input_records": 4,
            "surviving_records": 3,
            "rejected_records": 1,
            "rejected_topologies": {"2s2f_1|2s2f_1": 1},
            "derivative_sectors": {"1": 4},
            "completion_digest": "physical-completions",
            "source_sha256": "structural",
            "destination_sha256": "physical",
        },
        "democratic_models": 3,
        "species_models": 3,
        "propagator_models": 3,
        "historical": {
            "records": 3,
            "classes": 3,
            "reproduced": 3,
            "missing": 0,
            "exact_comparisons": 3 + seed,
            "artifact": {"path": "/legacy/op_1.dat", "sha256": "legacy"},
        },
        "topologies": {
            "generator": {"2s2f_1|2s2f_1": 8},
            "exact_classes": {"2s2f_1|2s2f_1": 3},
        },
        "completion_digests": {
            "generator": "raw",
            "structural_exact_classes": "structural",
            "exact_classes": "exact",
        },
        "round_trip_digests": {
            "generator": "raw-round-trip",
            "exact_classes": "exact-round-trip",
        },
        "artifacts": artifacts,
    }


def test_completed_report_validation_checks_every_artifact(tmp_path):
    report = sample_report(tmp_path / "operator", 0)
    report_path = tmp_path / "report.json"
    report_path.write_text(json.dumps(report), encoding="utf-8")
    inventory = {
        "operator": "1",
        "kind": "regular",
        "derivatives": 0,
        "legacy_sha256": "legacy",
    }

    assert completed_report_is_valid(report_path, inventory, 0)
    next(iter(report["artifacts"].values()))["sha256"] = "changed"
    report_path.write_text(json.dumps(report), encoding="utf-8")
    assert not completed_report_is_valid(report_path, inventory, 0)


def test_stable_report_comparison_ignores_only_operational_fields(tmp_path):
    seed_zero = sample_report(tmp_path / "seed-zero", 0)
    seed_one = sample_report(tmp_path / "seed-one", 1)

    assert stable_report_data(seed_zero) == stable_report_data(seed_one)
    seed_one["artifacts"]["generator"]["sha256"] = "order-dependent-raw"
    seed_one["artifacts"]["exact_classes"]["sha256"] = (
        "order-dependent-exact"
    )
    assert stable_report_data(seed_zero) == stable_report_data(seed_one)
    seed_one["completion_digests"]["generator"] = "changed"
    assert stable_report_data(seed_zero) != stable_report_data(seed_one)


def test_partial_manifest_requires_two_identical_seed_reports(tmp_path, monkeypatch):
    output_root = tmp_path / "rebuild"
    inventory = (
        {
            "operator": "1",
            "kind": "regular",
            "derivatives": 0,
            "legacy_path": "/legacy/op_1.dat",
            "legacy_sha256": "legacy",
        },
    )
    monkeypatch.setattr(
        "rebuild_completion_database.operator_inventory",
        lambda legacy_dir: inventory,
    )
    for seed in (0, 1):
        run_path = output_root / f"seed-{seed}" / "run.json"
        run_path.parent.mkdir(parents=True)
        run_path.write_text(
            json.dumps(
                {
                    "source_commit": "abc",
                    "hash_seed": str(seed),
                    "dirty_worktree": False,
                    "inventory": list(inventory),
                }
            ),
            encoding="utf-8",
        )
        report = sample_report(
            output_root / f"seed-{seed}" / "operators" / "1" / "artifacts",
            seed,
        )
        report_path = (
            output_root / f"seed-{seed}" / "operators" / "1" / "report.json"
        )
        report_path.write_text(json.dumps(report), encoding="utf-8")

    manifest = build_manifest(output_root, tmp_path / "legacy", operators=["1"])

    assert manifest["source_commit"] == "abc"
    assert manifest["manifest_operator_count"] == 1
    assert manifest["complete"]
    assert manifest["reproducible"]
    assert manifest["operators"]["1"]["census"]["records"]["generator"] == 8

    seed_one_path = output_root / "seed-1" / "operators" / "1" / "report.json"
    seed_one = json.loads(seed_one_path.read_text(encoding="utf-8"))
    seed_one["records"]["generator"] = 9
    seed_one_path.write_text(json.dumps(seed_one), encoding="utf-8")
    with pytest.raises(ValueError, match="Hash-seed mismatch for 1"):
        build_manifest(output_root, tmp_path / "legacy", operators=["1"])
