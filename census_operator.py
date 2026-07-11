#!/usr/bin/env python3

"""Regenerate and verify one operator in the full completion database."""

import argparse
import json
import os
from pathlib import Path
import tempfile

from neutrinomass.database.rebuild import census_operator, operator_registry


def write_report(path, report):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
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
            temporary.write(payload)
        temporary_path.replace(path)
    finally:
        temporary_path.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operator", choices=sorted(operator_registry()))
    parser.add_argument("historical_path", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    report = census_operator(
        args.operator,
        args.historical_path,
        args.output_dir,
        hash_seed=os.environ.get("PYTHONHASHSEED", ""),
    )
    write_report(args.report, report)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
