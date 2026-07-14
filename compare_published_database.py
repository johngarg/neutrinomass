#!/usr/bin/env python3

"""Compare a completed rebuild with the published 2009.13537 database."""

import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
import sqlite3
import tempfile

from neutrinomass.completions.fingerprints import lagrangian_fingerprint
from neutrinomass.database.rebuild import model_strings, read_model_artifact
from neutrinomass.database.serialization import loads_completion
from rebuild_completion_database import atomic_write_json, load_json


PUBLISHED_HEADLINE = {
    "physical_lagrangians": 430810,
    "democratic_models": 141989,
    "filtered_lagrangians": 11483,
    "filtered_democratic_models": 11216,
}
PUBLISHED_ONLY_OPERATORS = ("1ppp", "D1", "D11", "D19a", "D19b", "D19c")
CURRENT_ONLY_OPERATORS = ("D13c",)


def published_filtered_counts():
    """Return the per-operator counts in the packaged published database."""

    from neutrinomass.database.database import DATA

    return Counter(DATA["op"])


def read_filtered_models(path):
    result = {}
    with Path(path).open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, start=1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
                operator = record["operator"]
                fields = tuple(record["fields"])
            except (KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
                raise ValueError(
                    f"Invalid filtered-model record at {path}:{line_number}"
                ) from error
            result.setdefault(operator, set()).add(fields)
    return result


def build_operator_rows(manifest, filtering_report, published_filtered):
    rows = []
    current_only = set(CURRENT_ONLY_OPERATORS)
    for operator, entry in sorted(manifest["operators"].items()):
        census = entry["census"]
        filtering = filtering_report["operators"][operator]
        has_published_table_baseline = operator not in current_only
        published_models = (
            census["historical"]["records"]
            if has_published_table_baseline
            else None
        )
        published_survivors = (
            published_filtered.get(operator, 0)
            if has_published_table_baseline
            else None
        )
        updated_models = census["democratic_models"]
        updated_survivors = filtering["survivors"]
        rows.append(
            {
                "operator": operator,
                "kind": entry["kind"],
                "derivatives": entry["derivatives"],
                "legacy_records": census["historical"]["records"],
                "published_models": published_models,
                "updated_structural_classes": census[
                    "structural_exact_classes"
                ],
                "amplitude_rejected_classes": census[
                    "amplitude_symmetrisation"
                ]["rejected_records"],
                "updated_physical_classes": census["exact_classes"]["all"],
                "updated_species_models": census["species_models"],
                "updated_democratic_models": updated_models,
                "model_delta": (
                    updated_models - published_models
                    if published_models is not None
                    else None
                ),
                "published_filtered_models": published_survivors,
                "updated_filtered_models": updated_survivors,
                "filtered_delta": (
                    updated_survivors - published_survivors
                    if published_survivors is not None
                    else None
                ),
                "decoded_vertex_rejections": census[
                    "generation_rejections"
                ]["decoded_vanishing_uv_interactions"],
            }
        )
    return rows


def _iter_artifact_records(path):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, start=1):
            if not line.strip():
                continue
            try:
                yield loads_completion(line)
            except (AttributeError, IndexError, KeyError, TypeError, ValueError) as error:
                raise ValueError(
                    f"Invalid exact completion at {path}:{line_number}"
                ) from error


def scan_current_totals(output_root, manifest, filtered_models, work_dir):
    """Count global Lagrangians and models without retaining completions."""

    output_root = Path(output_root).resolve()
    input_models = set()
    filtered_model_union = set().union(*filtered_models.values())
    physical_incidences = 0
    filtered_incidences = 0
    database = tempfile.NamedTemporaryFile(
        prefix="published-comparison-",
        suffix=".sqlite3",
        dir=work_dir,
        delete=False,
    )
    database_path = Path(database.name)
    database.close()
    try:
        with sqlite3.connect(database_path) as connection:
            connection.execute(
                "CREATE TABLE physical (fingerprint TEXT PRIMARY KEY)"
            )
            connection.execute(
                "CREATE TABLE filtered (fingerprint TEXT PRIMARY KEY)"
            )
            for operator, entry in sorted(manifest["operators"].items()):
                model_path = output_root / entry["artifacts"][
                    "democratic_models"
                ]["path"]
                input_models.update(fields for fields, _ in read_model_artifact(model_path))
                exact_path = output_root / entry["artifacts"]["exact_classes"][
                    "path"
                ]
                operator_survivors = filtered_models.get(operator, set())
                for completion in _iter_artifact_records(exact_path):
                    physical_incidences += 1
                    fingerprint = repr(lagrangian_fingerprint(completion))
                    connection.execute(
                        "INSERT OR IGNORE INTO physical(fingerprint) VALUES (?)",
                        (fingerprint,),
                    )
                    if model_strings(completion) in operator_survivors:
                        filtered_incidences += 1
                        connection.execute(
                            "INSERT OR IGNORE INTO filtered(fingerprint) VALUES (?)",
                            (fingerprint,),
                        )
            connection.commit()
            physical = connection.execute(
                "SELECT COUNT(*) FROM physical"
            ).fetchone()[0]
            filtered = connection.execute(
                "SELECT COUNT(*) FROM filtered"
            ).fetchone()[0]
    finally:
        database_path.unlink(missing_ok=True)
    return {
        "physical_lagrangians": physical,
        "physical_lagrangian_operator_incidences": physical_incidences,
        "democratic_models": len(input_models),
        "filtered_lagrangians": filtered,
        "filtered_lagrangian_operator_incidences": filtered_incidences,
        "filtered_democratic_models": len(filtered_model_union),
    }


def _sum_rows(rows, key):
    return sum(row[key] for row in rows if row[key] is not None)


def build_report(manifest, filtering_report, rows, current_totals):
    comparable = [row for row in rows if row["published_models"] is not None]
    current = dict(current_totals)
    return {
        "schema_version": 1,
        "source_commit": manifest["source_commit"],
        "published_reference": {
            "paper": "2009.13537",
            "headline": PUBLISHED_HEADLINE,
            "published_only_operators": list(PUBLISHED_ONLY_OPERATORS),
            "current_only_operators": list(CURRENT_ONLY_OPERATORS),
        },
        "current": current,
        "headline_delta": {
            key: current[key] - value for key, value in PUBLISHED_HEADLINE.items()
        },
        "operator_incidences": {
            "comparable_operator_count": len(comparable),
            "published_models": _sum_rows(comparable, "published_models"),
            "updated_democratic_models": _sum_rows(
                comparable, "updated_democratic_models"
            ),
            "model_delta": _sum_rows(comparable, "model_delta"),
            "published_filtered_models": _sum_rows(
                comparable, "published_filtered_models"
            ),
            "updated_filtered_models": _sum_rows(
                comparable, "updated_filtered_models"
            ),
            "filtered_delta": _sum_rows(comparable, "filtered_delta"),
        },
        "filtering_totals": filtering_report["totals"],
        "operators": {row["operator"]: row for row in rows},
    }


def write_csv(path, rows):
    path = Path(path)
    fieldnames = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path, report, rows):
    labels = {
        "physical_lagrangians": "Inequivalent Lagrangians",
        "democratic_models": "Distinct unfiltered democratic models",
        "filtered_lagrangians": "Filtered Lagrangians",
        "filtered_democratic_models": "Distinct filtered democratic models",
    }
    lines = [
        "# Updated completion-database comparison",
        "",
        f"Scientific source commit: `{report['source_commit']}`",
        "",
        "| Quantity | Published | Updated | Delta |",
        "|---|---:|---:|---:|",
    ]
    for key, label in labels.items():
        lines.append(
            f"| {label} | {report['published_reference']['headline'][key]:,} "
            f"| {report['current'][key]:,} | {report['headline_delta'][key]:+,} |"
        )
    changed = [
        row
        for row in rows
        if row["model_delta"] not in (None, 0)
        or row["filtered_delta"] not in (None, 0)
    ]
    lines.extend(
        [
            "",
            "## Changed operator-level model counts",
            "",
            "| Operator | Published models | Updated models | Delta | "
            "Published filtered | Updated filtered | Delta |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in changed:
        lines.append(
            f"| {row['operator']} | {row['published_models']} | "
            f"{row['updated_democratic_models']} | {row['model_delta']:+} | "
            f"{row['published_filtered_models']} | "
            f"{row['updated_filtered_models']} | {row['filtered_delta']:+} |"
        )
    lines.extend(
        [
            "",
            "The headline comparison uses the four global counts quoted in "
            "2009.13537. The operator table comparison covers the common registry "
            "only. The published-only operators are "
            + ", ".join(PUBLISHED_ONLY_OPERATORS)
            + "; D13c is current-only.",
            "",
        ]
    )
    Path(path).write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_root", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("filtering_report", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = load_json(args.manifest)
    filtering_report = load_json(args.filtering_report)
    filtered_path = Path(filtering_report["filtered_models"]["path"])
    filtered_models = read_filtered_models(filtered_path)
    rows = build_operator_rows(
        manifest, filtering_report, published_filtered_counts()
    )
    current_totals = scan_current_totals(
        args.output_root, manifest, filtered_models, args.output_dir
    )
    report = build_report(manifest, filtering_report, rows, current_totals)
    json_path = args.output_dir / "published_comparison.json"
    csv_path = args.output_dir / "published_comparison.csv"
    markdown_path = args.output_dir / "published_comparison.md"
    atomic_write_json(json_path, report)
    write_csv(csv_path, rows)
    write_markdown(markdown_path, report, rows)
    print(json.dumps(report["headline_delta"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
