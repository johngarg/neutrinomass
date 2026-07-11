import gzip
import json
from pathlib import Path

from cluster_rebuild import (
    build_task_manifest,
    package_report,
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
