#!/usr/bin/env python3

"""Plan, run, and consolidate race-free Priority-4 cluster tasks."""

import argparse
from copy import deepcopy
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from time import time

from neutrinomass.database.rebuild import file_sha256, operator_inventory
from rebuild_completion_database import (
    REPOSITORY,
    atomic_write_json,
    build_manifest,
    completed_report_is_valid,
    load_json,
    operator_report_path,
    seed_directory,
)


SCHEMA_VERSION = 1
TASK_SCHEMA_VERSION = 2
SCIENTIFIC_PATHS = ("census_operator.py", "neutrinomass")
HIGH_RESOURCE_OPERATORS = frozenset(
    {
        "71p",
        "77p",
        "78p",
        "79a",
        "79b",
        "7p",
        "80a",
        "80b",
        "80c",
        "80d",
        "81a",
        "81b",
        "81c",
        "81d",
        "8pp",
    }
)


def git_output(*args):
    return subprocess.run(
        ["git", *args],
        cwd=REPOSITORY,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def validate_source_checkout(source_commit, orchestration_commit):
    """Require the planned orchestration and unchanged scientific sources."""

    if git_output("rev-parse", "HEAD") != orchestration_commit:
        raise ValueError("Cluster task manifest targets a different checkout")
    worktree = subprocess.run(
        ["git", "diff", "--quiet", "HEAD", "--"], cwd=REPOSITORY
    )
    if worktree.returncode == 1:
        raise ValueError("Cluster task checkout has modified tracked files")
    if worktree.returncode > 1:
        raise RuntimeError("Could not validate the cluster task checkout")
    result = subprocess.run(
        ["git", "diff", "--quiet", source_commit, "--", *SCIENTIFIC_PATHS],
        cwd=REPOSITORY,
    )
    if result.returncode == 1:
        raise ValueError(
            "Scientific census sources differ from the pinned source commit"
        )
    if result.returncode > 1:
        raise RuntimeError("Could not compare the scientific census sources")


def portable_inventory(legacy_dir):
    return [
        {
            "operator": item["operator"],
            "kind": item["kind"],
            "derivatives": item["derivatives"],
            "legacy_file": Path(item["legacy_path"]).name,
            "legacy_sha256": item["legacy_sha256"],
        }
        for item in operator_inventory(legacy_dir)
    ]


def operator_resource_tier(operator):
    return "high" if operator in HIGH_RESOURCE_OPERATORS else "normal"


def build_task_manifest(legacy_dir, source_commit, orchestration_commit):
    inventory = portable_inventory(legacy_dir)
    tasks = []
    for seed in (0, 1):
        for item in inventory:
            tasks.append(
                {
                    "task_id": len(tasks),
                    "hash_seed": str(seed),
                    "operator": item["operator"],
                    "resource_tier": operator_resource_tier(item["operator"]),
                }
            )
    return {
        "schema_version": TASK_SCHEMA_VERSION,
        "source_commit": source_commit,
        "orchestration_commit": orchestration_commit,
        "inventory": inventory,
        "tasks": tasks,
    }


def validate_task_manifest(manifest, legacy_dir, *, verify_inventory=True):
    if manifest.get("schema_version") != TASK_SCHEMA_VERSION:
        raise ValueError("Unsupported cluster task manifest schema")
    if verify_inventory and manifest["inventory"] != portable_inventory(legacy_dir):
        raise ValueError("Legacy inventory differs from the task manifest")
    expected = []
    for seed in (0, 1):
        for item in manifest["inventory"]:
            expected.append(
                {
                    "task_id": len(expected),
                    "hash_seed": str(seed),
                    "operator": item["operator"],
                    "resource_tier": operator_resource_tier(
                        item["operator"]
                    ),
                }
            )
    if manifest["tasks"] != expected:
        raise ValueError("Cluster task list is incomplete or out of order")


def task_array_expression(manifest, resource_tier):
    """Return a compact Slurm array expression for one resource tier."""

    if resource_tier not in {"normal", "high"}:
        raise ValueError(f"Unknown resource tier {resource_tier!r}")
    task_ids = [
        task["task_id"]
        for task in manifest["tasks"]
        if task["resource_tier"] == resource_tier
    ]
    if not task_ids:
        raise ValueError(f"No {resource_tier} resource tasks in manifest")

    ranges = []
    start = previous = task_ids[0]
    for task_id in task_ids[1:]:
        if task_id == previous + 1:
            previous = task_id
            continue
        ranges.append(str(start) if start == previous else f"{start}-{previous}")
        start = previous = task_id
    ranges.append(str(start) if start == previous else f"{start}-{previous}")
    return ",".join(ranges)


def _compressed_artifact(source, metadata):
    source = Path(source)
    if file_sha256(source) != metadata["sha256"]:
        raise ValueError(f"Artifact changed before compression: {source}")
    uncompressed_size = source.stat().st_size
    destination = source.with_suffix(source.suffix + ".gz")
    temporary = tempfile.NamedTemporaryFile(
        mode="w+b",
        prefix=f".{destination.name}.",
        suffix=".tmp",
        dir=source.parent,
        delete=False,
    )
    temporary_path = Path(temporary.name)
    try:
        with temporary:
            with gzip.GzipFile(
                filename="", mode="wb", fileobj=temporary, mtime=0
            ) as compressed:
                with source.open("rb") as stream:
                    shutil.copyfileobj(stream, compressed, length=1024 * 1024)
        compressed_sha256 = file_sha256(temporary_path)
        digest = hashlib.sha256()
        size = 0
        with gzip.open(temporary_path, "rb") as stream:
            while True:
                block = stream.read(1024 * 1024)
                if not block:
                    break
                digest.update(block)
                size += len(block)
        if digest.hexdigest() != metadata["sha256"] or size != uncompressed_size:
            raise ValueError(f"Compression round trip failed: {source}")
        temporary_path.replace(destination)
    finally:
        temporary_path.unlink(missing_ok=True)
    return destination, {
        "path": str(destination.resolve()),
        "sha256": compressed_sha256,
        "compression": "gzip",
        "uncompressed_sha256": metadata["sha256"],
        "uncompressed_size_bytes": uncompressed_size,
        "compressed_size_bytes": destination.stat().st_size,
    }


def package_report(report_path):
    """Compress large safe artifacts and update the report atomically."""

    report_path = Path(report_path)
    report = load_json(report_path)
    originals = []
    for name in (
        "generator",
        "structural_exact_classes",
        "exact_classes",
    ):
        metadata = report["artifacts"][name]
        source = Path(metadata["path"])
        destination, compressed = _compressed_artifact(source, metadata)
        originals.append(source)
        report["artifacts"][name] = compressed
        if not destination.exists():
            raise ValueError(f"Missing compressed artifact: {destination}")
    atomic_write_json(report_path, report)
    for source in originals:
        source.unlink()
    return report


def _stage_report(report, scratch_dir, final_dir, provenance):
    """Copy a packaged task to shared storage, publishing report.json last."""

    scratch_dir = Path(scratch_dir)
    final_dir = Path(final_dir).resolve()
    final_dir.parent.mkdir(parents=True, exist_ok=True)
    staging_dir = Path(
        tempfile.mkdtemp(prefix=f".{final_dir.name}.stage-", dir=final_dir.parent)
    )
    staged_report = deepcopy(report)
    try:
        for name, metadata in staged_report["artifacts"].items():
            source = Path(metadata["path"])
            destination = staging_dir / source.name
            shutil.copy2(source, destination)
            if file_sha256(destination) != metadata["sha256"]:
                raise ValueError(f"Staged artifact checksum failed: {name}")
            metadata["path"] = str((final_dir / source.name).resolve())

        log_path = scratch_dir / "census.log"
        if log_path.exists():
            shutil.copy2(log_path, staging_dir / log_path.name)

        raw_name = f"op_{staged_report['operator']}_remediated.jsonl"
        exact_name = (
            f"op_{staged_report['operator']}_remediated_unique.jsonl"
        )
        staged_report["deduplication"]["source"] = str(final_dir / raw_name)
        staged_report["deduplication"]["destination"] = str(
            final_dir / exact_name
        )
        staged_report["cluster_task"] = provenance
        atomic_write_json(staging_dir / "report.json", staged_report)

        final_dir.mkdir(parents=True, exist_ok=True)
        for path in staging_dir.iterdir():
            if path.name != "report.json":
                path.replace(final_dir / path.name)
        (staging_dir / "report.json").replace(final_dir / "report.json")
    finally:
        shutil.rmtree(staging_dir, ignore_errors=True)
    return staged_report


def _report_source_is_valid(report_path, output_root, seed, source_commit):
    report = load_json(report_path)
    provenance = report.get("cluster_task")
    if provenance is not None:
        return provenance.get("source_commit") == source_commit
    run_path = seed_directory(output_root, seed) / "run.json"
    return (
        run_path.exists()
        and load_json(run_path).get("source_commit") == source_commit
    )


def relocate_report_paths(report_path, legacy_path):
    """Make a transferred self-contained report refer to its new directory."""

    report_path = Path(report_path)
    if not report_path.exists():
        return False
    report = load_json(report_path)
    changed = False
    for artifact in report["artifacts"].values():
        recorded = Path(artifact["path"])
        candidate = (report_path.parent / recorded.name).resolve()
        if recorded != candidate and candidate.exists():
            if file_sha256(candidate) != artifact["sha256"]:
                raise ValueError(f"Relocated artifact checksum failed: {candidate}")
            artifact["path"] = str(candidate)
            changed = True

    legacy_path = Path(legacy_path).resolve()
    historical = report["historical"]["artifact"]
    if Path(historical["path"]) != legacy_path:
        historical["path"] = str(legacy_path)
        changed = True

    operator = report["operator"]
    raw_path = report_path.parent / f"op_{operator}_remediated.jsonl"
    exact_path = report_path.parent / f"op_{operator}_remediated_unique.jsonl"
    if Path(report["deduplication"]["source"]) != raw_path:
        report["deduplication"]["source"] = str(raw_path.resolve())
        changed = True
    if Path(report["deduplication"]["destination"]) != exact_path:
        report["deduplication"]["destination"] = str(exact_path.resolve())
        changed = True
    if changed:
        atomic_write_json(report_path, report)
    return changed


def run_task(manifest_path, task_id, output_root, legacy_dir, scratch_root):
    manifest = load_json(manifest_path)
    # The plan and final consolidation hash all 243 legacy inputs. A worker
    # hashes only its own input so 486 array jobs do not reread the archive.
    validate_task_manifest(manifest, legacy_dir, verify_inventory=False)
    validate_source_checkout(
        manifest["source_commit"], manifest["orchestration_commit"]
    )
    try:
        task = manifest["tasks"][task_id]
    except IndexError as error:
        raise ValueError(f"Unknown cluster task {task_id}") from error
    if task["task_id"] != task_id:
        raise ValueError("Task index does not match its manifest entry")

    by_name = {item["operator"]: item for item in manifest["inventory"]}
    item = by_name[task["operator"]]
    inventory_item = {
        "operator": item["operator"],
        "kind": item["kind"],
        "derivatives": item["derivatives"],
        "legacy_sha256": item["legacy_sha256"],
    }
    report_path = operator_report_path(
        output_root, task["hash_seed"], task["operator"]
    )
    legacy_path = Path(legacy_dir) / item["legacy_file"]
    if file_sha256(legacy_path) != item["legacy_sha256"]:
        raise ValueError(f"Legacy checksum mismatch: {legacy_path}")
    relocate_report_paths(report_path, legacy_path)
    if completed_report_is_valid(
        report_path, inventory_item, task["hash_seed"]
    ) and _report_source_is_valid(
        report_path,
        output_root,
        task["hash_seed"],
        manifest["source_commit"],
    ):
        return {"status": "already_complete", "report": str(report_path)}

    scratch_root = Path(scratch_root)
    scratch_root.mkdir(parents=True, exist_ok=True)
    started = time()
    with tempfile.TemporaryDirectory(
        prefix=f"priority4-{task['hash_seed']}-{task['operator']}-",
        dir=scratch_root,
    ) as temporary:
        scratch_dir = Path(temporary)
        scratch_report = scratch_dir / "report.json"
        environment = os.environ.copy()
        environment["PYTHONHASHSEED"] = task["hash_seed"]
        environment["PYTHONWARNINGS"] = "ignore::SyntaxWarning"
        environment["MPLCONFIGDIR"] = str(scratch_dir / ".matplotlib")
        command = [
            sys.executable,
            str(REPOSITORY / "census_operator.py"),
            task["operator"],
            str(legacy_path),
            str(scratch_dir),
            "--report",
            str(scratch_report),
        ]
        with (scratch_dir / "census.log").open("w", encoding="utf-8") as log:
            result = subprocess.run(
                command,
                cwd=REPOSITORY,
                env=environment,
                stdout=log,
                stderr=subprocess.STDOUT,
            )
        if result.returncode:
            print(
                (scratch_dir / "census.log").read_text(
                    encoding="utf-8", errors="replace"
                ),
                file=sys.stderr,
            )
            raise RuntimeError(
                f"Task {task_id} failed with exit code {result.returncode}"
            )
        report = package_report(scratch_report)
        provenance = {
            "schema_version": SCHEMA_VERSION,
            "task_id": task_id,
            "source_commit": manifest["source_commit"],
            "orchestration_commit": manifest["orchestration_commit"],
            "started_unix": started,
            "finished_unix": time(),
        }
        final_dir = report_path.parent
        _stage_report(report, scratch_dir, final_dir, provenance)

    if not completed_report_is_valid(
        report_path, inventory_item, task["hash_seed"]
    ):
        raise ValueError(f"Task {task_id} staged an invalid report")
    return {"status": "complete", "report": str(report_path)}


def consolidate(output_root, legacy_dir, task_manifest, manifest_output):
    task_data = load_json(task_manifest)
    validate_task_manifest(task_data, legacy_dir)
    complete_inventory = operator_inventory(legacy_dir)
    by_name = {item["operator"]: item for item in complete_inventory}

    for seed in (0, 1):
        original_run_path = seed_directory(output_root, seed) / "run.json"
        original_run = (
            load_json(original_run_path) if original_run_path.exists() else None
        )
        operators = {}
        for item in complete_inventory:
            name = item["operator"]
            report_path = operator_report_path(output_root, seed, name)
            if not completed_report_is_valid(report_path, item, seed):
                raise ValueError(f"Invalid or missing seed-{seed} report for {name}")
            report = load_json(report_path)
            provenance = report.get("cluster_task")
            cluster_source_matches = provenance is not None and (
                provenance.get("source_commit") == task_data["source_commit"]
            )
            serial_source_matches = provenance is None and original_run and (
                original_run.get("source_commit") == task_data["source_commit"]
            )
            if not (cluster_source_matches or serial_source_matches):
                raise ValueError(f"Source-commit mismatch for seed-{seed} {name}")
            operators[name] = {
                "status": "complete",
                "report": str(report_path.resolve()),
            }
        run = {
            "schema_version": SCHEMA_VERSION,
            "hash_seed": str(seed),
            "source_commit": task_data["source_commit"],
            "orchestration_commit": task_data["orchestration_commit"],
            "dirty_worktree": False,
            "inventory": list(complete_inventory),
            "selected_operators": list(by_name),
            "operators": operators,
        }
        atomic_write_json(seed_directory(output_root, seed) / "run.json", run)

    migration = build_manifest(output_root, legacy_dir)
    atomic_write_json(manifest_output, migration)
    return migration


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    plan = subparsers.add_parser("plan", help="write the 486-task manifest")
    plan.add_argument("legacy_dir", type=Path)
    plan.add_argument("--source-commit", required=True)
    plan.add_argument("--output", type=Path, required=True)

    count = subparsers.add_parser("count", help="print the task count")
    count.add_argument("task_manifest", type=Path)

    array = subparsers.add_parser(
        "array", help="print the Slurm array expression for a resource tier"
    )
    array.add_argument("task_manifest", type=Path)
    array.add_argument("resource_tier", choices=("normal", "high"))

    worker = subparsers.add_parser("worker", help="run one immutable task")
    worker.add_argument("task_manifest", type=Path)
    worker.add_argument("task_id", type=int)
    worker.add_argument("output_root", type=Path)
    worker.add_argument("legacy_dir", type=Path)
    worker.add_argument("--scratch-root", type=Path, required=True)

    merge = subparsers.add_parser(
        "consolidate", help="validate all tasks and create seed/run manifests"
    )
    merge.add_argument("task_manifest", type=Path)
    merge.add_argument("output_root", type=Path)
    merge.add_argument("legacy_dir", type=Path)
    merge.add_argument("--output", type=Path, required=True)

    args = parser.parse_args()
    if args.command == "plan":
        source_commit = git_output("rev-parse", args.source_commit)
        orchestration_commit = git_output("rev-parse", "HEAD")
        validate_source_checkout(source_commit, orchestration_commit)
        manifest = build_task_manifest(
            args.legacy_dir, source_commit, orchestration_commit
        )
        atomic_write_json(args.output, manifest)
        print(json.dumps({"tasks": len(manifest["tasks"])}, sort_keys=True))
        return
    if args.command == "count":
        print(len(load_json(args.task_manifest)["tasks"]))
        return
    if args.command == "array":
        print(
            task_array_expression(
                load_json(args.task_manifest), args.resource_tier
            )
        )
        return
    if args.command == "worker":
        result = run_task(
            args.task_manifest,
            args.task_id,
            args.output_root,
            args.legacy_dir,
            args.scratch_root,
        )
        print(json.dumps(result, indent=2, sort_keys=True))
        return
    result = consolidate(
        args.output_root,
        args.legacy_dir,
        args.task_manifest,
        args.output,
    )
    print(
        json.dumps(
            {
                "complete": result["complete"],
                "operators": result["manifest_operator_count"],
                "reproducible": result["reproducible"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
