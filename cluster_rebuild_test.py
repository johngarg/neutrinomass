import gzip
import json
from pathlib import Path

from cluster_rebuild import (
    build_task_manifest,
    package_report,
    relocate_report_paths,
    validate_task_manifest,
)
from neutrinomass.database.rebuild import file_sha256


def test_task_manifest_is_seed_major_and_complete(tmp_path, monkeypatch):
    inventory = [
        {
            "operator": "1",
            "kind": "regular",
            "derivatives": 0,
            "legacy_path": str(tmp_path / "op_1.dat"),
            "legacy_sha256": "one",
        },
        {
            "operator": "1p",
            "kind": "regular",
            "derivatives": 0,
            "legacy_path": str(tmp_path / "op_1p.dat"),
            "legacy_sha256": "one-prime",
        },
    ]
    monkeypatch.setattr("cluster_rebuild.operator_inventory", lambda _: inventory)

    manifest = build_task_manifest(tmp_path, "science", "cluster")

    assert manifest["source_commit"] == "science"
    assert [(task["hash_seed"], task["operator"]) for task in manifest["tasks"]] == [
        ("0", "1"),
        ("0", "1p"),
        ("1", "1"),
        ("1", "1p"),
    ]
    validate_task_manifest(manifest, tmp_path)


def test_worker_manifest_validation_does_not_rehash_entire_legacy_archive(
    tmp_path, monkeypatch
):
    manifest = {
        "schema_version": 1,
        "inventory": [
            {
                "operator": "1",
                "kind": "regular",
                "derivatives": 0,
                "legacy_file": "op_1.dat",
                "legacy_sha256": "one",
            }
        ],
        "tasks": [{"task_id": 0, "hash_seed": "0", "operator": "1"},
                  {"task_id": 1, "hash_seed": "1", "operator": "1"}],
    }
    monkeypatch.setattr(
        "cluster_rebuild.portable_inventory",
        lambda _: (_ for _ in ()).throw(AssertionError("rehash")),
    )

    validate_task_manifest(manifest, tmp_path, verify_inventory=False)


def test_package_report_is_deterministic_and_lossless(tmp_path):
    artifacts = {}
    original_payloads = {}
    for name in ("generator", "exact_classes", "democratic_models"):
        path = tmp_path / f"{name}.jsonl"
        payload = ((name + "\n") * 100).encode()
        path.write_bytes(payload)
        original_payloads[name] = payload
        artifacts[name] = {"path": str(path), "sha256": file_sha256(path)}
    report_path = tmp_path / "report.json"
    report_path.write_text(json.dumps({"artifacts": artifacts}), encoding="utf-8")

    report = package_report(report_path)

    for name in ("generator", "exact_classes"):
        artifact = report["artifacts"][name]
        assert artifact["compression"] == "gzip"
        assert artifact["uncompressed_sha256"] == artifacts[name]["sha256"]
        assert not Path(artifacts[name]["path"]).exists()
        with gzip.open(artifact["path"], "rb") as stream:
            assert stream.read() == original_payloads[name]
    assert report["artifacts"]["democratic_models"] == artifacts[
        "democratic_models"
    ]

    first_hashes = {
        name: report["artifacts"][name]["sha256"]
        for name in ("generator", "exact_classes")
    }
    for name in ("generator", "exact_classes"):
        path = tmp_path / f"second-{name}.jsonl"
        path.write_bytes(original_payloads[name])
        artifacts[name] = {"path": str(path), "sha256": file_sha256(path)}
    second_report_path = tmp_path / "second-report.json"
    second_report_path.write_text(
        json.dumps({"artifacts": artifacts}), encoding="utf-8"
    )
    second = package_report(second_report_path)
    assert first_hashes == {
        name: second["artifacts"][name]["sha256"]
        for name in ("generator", "exact_classes")
    }


def test_relocate_report_paths_verifies_transferred_artifacts(tmp_path):
    old_root = Path("/old/machine/rebuild/seed-0/operators/1")
    artifacts = {}
    for name in ("generator", "exact_classes", "democratic_models"):
        path = tmp_path / f"{name}.jsonl"
        path.write_text(name, encoding="utf-8")
        artifacts[name] = {
            "path": str(old_root / path.name),
            "sha256": file_sha256(path),
        }
    legacy = tmp_path / "op_1.dat"
    legacy.write_text("legacy", encoding="utf-8")
    report_path = tmp_path / "report.json"
    report_path.write_text(
        json.dumps(
            {
                "operator": "1",
                "artifacts": artifacts,
                "historical": {
                    "artifact": {
                        "path": "/old/machine/raw_completions/op_1.dat",
                        "sha256": file_sha256(legacy),
                    }
                },
                "deduplication": {
                    "source": str(old_root / "op_1_remediated.jsonl"),
                    "destination": str(
                        old_root / "op_1_remediated_unique.jsonl"
                    ),
                },
            }
        ),
        encoding="utf-8",
    )

    assert relocate_report_paths(report_path, legacy)
    relocated = json.loads(report_path.read_text(encoding="utf-8"))
    assert all(
        Path(artifact["path"]).parent == tmp_path
        for artifact in relocated["artifacts"].values()
    )
    assert relocated["historical"]["artifact"]["path"] == str(
        legacy.resolve()
    )
