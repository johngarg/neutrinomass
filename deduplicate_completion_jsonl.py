#!/usr/bin/env python3

"""Deduplicate a safe completion JSONL artifact with bounded Python memory."""

import argparse
import json
from pathlib import Path

from neutrinomass.database import deduplicate_completion_jsonl


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()

    report = deduplicate_completion_jsonl(
        args.source,
        args.destination,
        work_dir=args.work_dir,
    )
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(payload, encoding="utf-8")
    print(payload, end="")


if __name__ == "__main__":
    main()
