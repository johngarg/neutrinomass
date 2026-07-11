#!/usr/bin/env python3

"""Run and verify the restartable full completion-database rebuild."""

import argparse
from collections import Counter
from copy import deepcopy
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from time import time

from neutrinomass.database.rebuild import (
    file_sha256,
    filter_democratic_registry,
    generated_one_loop_weinberg_models,
    operator_inventory,
    read_model_artifact,
    write_filtered_registry,
)


REPOSITORY = Path(__file__).resolve().parent
CENSUS_SCRIPT = REPOSITORY / "census_operator.py"


def atomic_write_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
        delete=False,
    )
    temporary_path = Path(temporary.name)
    try:
        with temporary:
            json.dump(payload, temporary, indent=2, sort_keys=True)
            temporary.write("\n")
        temporary_path.replace(path)
    finally:
        temporary_path.unlink(missing_ok=True)


def load_json(path):
    with Path(path).open(encoding="utf-8") as stream:
        return json.load(stream)


def git_state(repository=REPOSITORY):
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repository,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    dirty = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=repository,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    return commit, bool(dirty)


def seed_directory(output_root, seed):
    return Path(output_root) / f"seed-{seed}"


def operator_directory(output_root, seed, operator_name):
    return seed_directory(output_root, seed) / "operators" / operator_name


def operator_report_path(output_root, seed, operator_name):
    return operator_directory(output_root, seed, operator_name) / "report.json"


def completed_report_is_valid(report_path, inventory_item, seed):
    report_path = Path(report_path)
    if not report_path.exists():
        return False
    try:
        report = load_json(report_path)
        if report["schema_version"] != 1:
            return False
        if report["operator"] != inventory_item["operator"]:
            return False
        if report["kind"] != inventory_item["kind"]:
            return False
        if report["derivatives"] != inventory_item["derivatives"]:
            return False
        if report["hash_seed"] != str(seed):
            return False
        historical = report["historical"]["artifact"]
        if historical["sha256"] != inventory_item["legacy_sha256"]:
            return False
        for artifact in report["artifacts"].values():
            path = Path(artifact["path"])
            if not path.exists() or file_sha256(path) != artifact["sha256"]:
                return False
    except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError):
        return False
    return True


def _initial_run_metadata(seed, commit, dirty, inventory, selected):
    return {
        "schema_version": 1,
        "hash_seed": str(seed),
        "source_commit": commit,
        "dirty_worktree": dirty,
        "inventory": list(inventory),
        "selected_operators": list(selected),
        "operators": {
            name: {"status": "pending"} for name in selected
        },
    }


def run_seed(
    output_root,
    legacy_dir,
    seed,
    *,
    operators=None,
    resume=False,
    allow_dirty=False,
):
    output_root = Path(output_root).resolve()
    inventory = operator_inventory(legacy_dir)
    by_name = {item["operator"]: item for item in inventory}
    selected = tuple(operators or by_name)
    unknown = sorted(set(selected) - set(by_name))
    if unknown:
        raise ValueError(f"Unknown operators: {unknown}")

    commit, dirty = git_state()
    if dirty and not allow_dirty:
        raise ValueError(
            "Refusing a reproducibility run from a dirty worktree; "
            "commit the rebuild code or pass --allow-dirty for a pilot"
        )

    run_path = seed_directory(output_root, seed) / "run.json"
    if resume and run_path.exists():
        metadata = load_json(run_path)
        if metadata["source_commit"] != commit:
            raise ValueError("Cannot resume a seed run from a different commit")
        if metadata["hash_seed"] != str(seed):
            raise ValueError("Seed metadata does not match the requested seed")
        if metadata["inventory"] != list(inventory):
            raise ValueError("Legacy inventory changed since the seed run began")
        missing = [name for name in selected if name not in metadata["operators"]]
        for name in missing:
            metadata["operators"][name] = {"status": "pending"}
        metadata["selected_operators"] = sorted(metadata["operators"])
    else:
        metadata = _initial_run_metadata(
            seed, commit, dirty, inventory, selected
        )
    atomic_write_json(run_path, metadata)

    for operator_name in selected:
        item = by_name[operator_name]
        report_path = operator_report_path(output_root, seed, operator_name)
        if resume and completed_report_is_valid(report_path, item, seed):
            metadata["operators"][operator_name] = {
                "status": "complete",
                "report": str(report_path),
                "resumed": True,
            }
            atomic_write_json(run_path, metadata)
            continue

        directory = operator_directory(output_root, seed, operator_name)
        directory.mkdir(parents=True, exist_ok=True)
        log_path = directory / "census.log"
        metadata["operators"][operator_name] = {
            "status": "running",
            "started_unix": time(),
            "report": str(report_path),
            "log": str(log_path),
        }
        atomic_write_json(run_path, metadata)
        command = [
            sys.executable,
            str(CENSUS_SCRIPT),
            operator_name,
            item["legacy_path"],
            str(directory),
            "--report",
            str(report_path),
        ]
        environment = os.environ.copy()
        environment["PYTHONHASHSEED"] = str(seed)
        environment["PYTHONWARNINGS"] = "ignore::SyntaxWarning"
        environment["MPLCONFIGDIR"] = str(directory / ".matplotlib")
        with log_path.open("w", encoding="utf-8") as log:
            result = subprocess.run(
                command,
                cwd=REPOSITORY,
                env=environment,
                stdout=log,
                stderr=subprocess.STDOUT,
            )
        if result.returncode:
            metadata["operators"][operator_name] = {
                **metadata["operators"][operator_name],
                "status": "failed",
                "returncode": result.returncode,
                "finished_unix": time(),
            }
            atomic_write_json(run_path, metadata)
            raise RuntimeError(
                f"{operator_name} failed with exit code {result.returncode}; "
                f"see {log_path}"
            )
        if not completed_report_is_valid(report_path, item, seed):
            raise ValueError(f"{operator_name} produced an invalid report")
        metadata["operators"][operator_name] = {
            **metadata["operators"][operator_name],
            "status": "complete",
            "finished_unix": time(),
        }
        atomic_write_json(run_path, metadata)
    return metadata


def stable_report_data(report):
    historical = deepcopy(report["historical"])
    historical.pop("exact_comparisons", None)
    historical["artifact"] = {
        "sha256": historical["artifact"]["sha256"]
    }
    return {
        "schema_version": report["schema_version"],
        "operator": report["operator"],
        "kind": report["kind"],
        "derivatives": report["derivatives"],
        "operator_scale_gev": report["operator_scale_gev"],
        "records": report["records"],
        "exact_classes": report["exact_classes"],
        "democratic_models": report["democratic_models"],
        "historical": historical,
        "topologies": report["topologies"],
        "completion_digests": report["completion_digests"],
        "round_trip_digests": report["round_trip_digests"],
        "democratic_model_sha256": report["artifacts"][
            "democratic_models"
        ]["sha256"],
    }


def build_manifest(output_root, legacy_dir, *, operators=None):
    output_root = Path(output_root).resolve()
    complete_inventory = operator_inventory(legacy_dir)
    by_name = {item["operator"]: item for item in complete_inventory}
    selected = tuple(operators or by_name)
    unknown = sorted(set(selected) - set(by_name))
    if unknown:
        raise ValueError(f"Unknown operators: {unknown}")
    inventory = tuple(by_name[name] for name in selected)
    seed_runs = {
        seed: load_json(seed_directory(output_root, seed) / "run.json")
        for seed in (0, 1)
    }
    commits = {run["source_commit"] for run in seed_runs.values()}
    if len(commits) != 1:
        raise ValueError("Seed runs used different source commits")
    for seed, run in seed_runs.items():
        if run["hash_seed"] != str(seed):
            raise ValueError(f"Run metadata has the wrong seed for seed {seed}")
        if run["inventory"] != list(complete_inventory):
            raise ValueError(f"Legacy inventory changed for seed {seed}")
    reproducible = not any(
        run.get("dirty_worktree", True) for run in seed_runs.values()
    )
    if len(inventory) == len(complete_inventory) and not reproducible:
        raise ValueError("Complete manifests require clean-worktree seed runs")

    entries = {}
    for item in inventory:
        operator_name = item["operator"]
        for seed in (0, 1):
            if not completed_report_is_valid(
                operator_report_path(output_root, seed, operator_name),
                item,
                seed,
            ):
                raise ValueError(
                    f"Invalid or corrupted seed-{seed} report for {operator_name}"
                )
        reports = {
            seed: load_json(
                operator_report_path(output_root, seed, operator_name)
            )
            for seed in (0, 1)
        }
        stable = {seed: stable_report_data(report) for seed, report in reports.items()}
        if stable[0] != stable[1]:
            raise ValueError(f"Hash-seed mismatch for {operator_name}")
        entries[operator_name] = {
            "kind": item["kind"],
            "derivatives": item["derivatives"],
            "legacy": {
                "path": str(
                    Path(Path(legacy_dir).name) / Path(item["legacy_path"]).name
                ),
                "sha256": item["legacy_sha256"],
            },
            "artifacts": {
                name: {
                    "path": str(
                        Path(artifact["path"]).resolve().relative_to(output_root)
                    ),
                    "sha256": artifact["sha256"],
                }
                for name, artifact in reports[0]["artifacts"].items()
            },
            "census": stable[0],
            "seed_reports": {
                str(seed): str(
                    operator_report_path(output_root, seed, operator_name)
                    .resolve()
                    .relative_to(output_root)
                )
                for seed in (0, 1)
            },
            "seed_artifact_sha256": {
                str(seed): {
                    name: artifact["sha256"]
                    for name, artifact in reports[seed]["artifacts"].items()
                }
                for seed in (0, 1)
            },
        }
    return {
        "schema_version": 1,
        "source_commit": commits.pop(),
        "hash_seeds": ["0", "1"],
        "legacy_operator_count": len(complete_inventory),
        "manifest_operator_count": len(inventory),
        "complete": len(inventory) == len(complete_inventory),
        "reproducible": reproducible,
        "operators": entries,
    }


def write_one_loop_models(path, models):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
        delete=False,
    )
    temporary_path = Path(temporary.name)
    try:
        with temporary:
            for fields in models:
                temporary.write(
                    json.dumps(
                        {"fields": list(fields)},
                        sort_keys=True,
                        separators=(",", ":"),
                    )
                    + "\n"
                )
        temporary_path.replace(path)
    finally:
        temporary_path.unlink(missing_ok=True)
    return file_sha256(path)


def filter_rebuilt_registry(output_root, manifest, output_dir):
    output_root = Path(output_root).resolve()
    output_dir = Path(output_dir).resolve()
    registry = {}
    scales = {}
    input_counts = {}
    for operator, entry in manifest["operators"].items():
        model_path = output_root / entry["artifacts"]["democratic_models"]["path"]
        models = dict(read_model_artifact(model_path))
        expected = entry["census"]["democratic_models"]
        if len(models) != expected:
            raise ValueError(
                f"Democratic-model count mismatch for {operator}: "
                f"expected {expected}, found {len(models)}"
            )
        registry[operator] = models
        scales[operator] = entry["census"]["operator_scale_gev"]
        input_counts[operator] = len(models)

    one_loop_models = generated_one_loop_weinberg_models()
    filtered = filter_democratic_registry(registry, scales, one_loop_models)
    survivor_path = output_dir / "filtered_democratic_models.jsonl"
    survivor_artifact = write_filtered_registry(
        survivor_path, filtered["survivors"]
    )
    one_loop_path = output_dir / "one_loop_weinberg_models.jsonl"
    one_loop_sha256 = write_one_loop_models(one_loop_path, one_loop_models)
    per_operator = {}
    for operator in sorted(registry):
        removed_mass = filtered["removed_by_mass"].get(operator, 0)
        removed_loop = filtered["removed_by_one_loop_weinberg"].get(operator, 0)
        survivors = len(filtered["survivors"][operator])
        if input_counts[operator] - removed_mass - removed_loop != survivors:
            raise ValueError(f"Filtering counts do not close for {operator}")
        per_operator[operator] = {
            "input": input_counts[operator],
            "removed_by_mass": removed_mass,
            "removed_by_one_loop_weinberg": removed_loop,
            "survivors": survivors,
        }
    return {
        "schema_version": 1,
        "source_commit": manifest["source_commit"],
        "complete_registry": manifest["complete"],
        "operator_count": len(registry),
        "one_loop_scale_gev": filtered["one_loop_scale_gev"],
        "one_loop_models": {
            "records": len(one_loop_models),
            "path": str(one_loop_path),
            "sha256": one_loop_sha256,
        },
        "filtered_models": {
            "records": survivor_artifact["records"],
            "path": str(survivor_path),
            "sha256": survivor_artifact["sha256"],
        },
        "totals": {
            "input": sum(input_counts.values()),
            "removed_by_mass": sum(filtered["removed_by_mass"].values()),
            "removed_by_one_loop_weinberg": sum(
                filtered["removed_by_one_loop_weinberg"].values()
            ),
            "survivors": survivor_artifact["records"],
        },
        "operators": per_operator,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="run one hash-seed rebuild")
    run_parser.add_argument("output_root", type=Path)
    run_parser.add_argument("legacy_dir", type=Path)
    run_parser.add_argument("--seed", choices=("0", "1"), required=True)
    run_parser.add_argument("--operators", nargs="*")
    run_parser.add_argument("--resume", action="store_true")
    run_parser.add_argument("--allow-dirty", action="store_true")

    manifest_parser = subparsers.add_parser(
        "manifest", help="compare seeds and write the migration manifest"
    )
    manifest_parser.add_argument("output_root", type=Path)
    manifest_parser.add_argument("legacy_dir", type=Path)
    manifest_parser.add_argument("--output", type=Path, required=True)
    manifest_parser.add_argument("--operators", nargs="*")

    filter_parser = subparsers.add_parser(
        "filter", help="filter the regenerated democratic model registry"
    )
    filter_parser.add_argument("output_root", type=Path)
    filter_parser.add_argument("manifest", type=Path)
    filter_parser.add_argument("--output-dir", type=Path, required=True)

    args = parser.parse_args()
    if args.command == "run":
        result = run_seed(
            args.output_root,
            args.legacy_dir,
            args.seed,
            operators=args.operators,
            resume=args.resume,
            allow_dirty=args.allow_dirty,
        )
        statuses = Counter(
            item["status"] for item in result["operators"].values()
        )
        print(
            json.dumps(
                {
                    "hash_seed": result["hash_seed"],
                    "source_commit": result["source_commit"],
                    "selected_operators": len(result["selected_operators"]),
                    "statuses": dict(sorted(statuses.items())),
                },
                indent=2,
                sort_keys=True,
            )
        )
        return

    if args.command == "manifest":
        manifest = build_manifest(
            args.output_root, args.legacy_dir, operators=args.operators
        )
        atomic_write_json(args.output, manifest)
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return

    manifest = load_json(args.manifest)
    report = filter_rebuilt_registry(
        args.output_root, manifest, args.output_dir
    )
    report_path = args.output_dir / "filtering_report.json"
    atomic_write_json(report_path, report)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
