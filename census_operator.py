#!/usr/bin/env python3

"""Regenerate and verify one operator in the full completion database."""

import argparse
import json
import os
from pathlib import Path
import tempfile
import sys

from neutrinomass.database.rebuild import census_operator, operator_registry
from neutrinomass.database.checkpoints import WorkBudget,BudgetExpired


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
    parser.add_argument("--report", type=Path)
    parser.add_argument("--stage",choices=("all","generation","decoded","local-reduction","structural","amplitudes","historical","coverage","package"),default="all")
    parser.add_argument("--unit",type=int,help="Stable zero-based generation-plan unit")
    parser.add_argument("--stage-unit",type=int,help="Stable 128-class amplitude or historical-coverage shard")
    parser.add_argument("--topology",help="Select matching topology units from the stable plan")
    parser.add_argument("--assignment-start",type=int,default=0)
    parser.add_argument("--assignment-stop",type=int)
    parser.add_argument("--producer",choices=("regular","local","routed"))
    parser.add_argument("--placement",type=int)
    parser.add_argument("--colour",type=int)
    parser.add_argument("--assignment-size",type=int,default=256)
    parser.add_argument("--resume",action="store_true")
    parser.add_argument("--time-budget",type=float,default=5400,help="Seconds before a saved continuation (exit 75)")
    parser.add_argument("--scratch",type=Path)
    parser.add_argument("--allow-dirty",action="store_true")
    parser.add_argument("--fault-after",type=int,help="Validation only: checkpoint and exit after N safe boundaries")
    parser.add_argument("--require-slurm",action="store_true")
    args = parser.parse_args()
    if args.require_slurm and not os.environ.get("SLURM_JOB_ID"):
        parser.error("Computations on Spartan require a Slurm allocation")
    budget=WorkBudget(args.time_budget,fault_after=args.fault_after)
    budget.scratch_root=args.scratch
    try:
        with budget.signals():
            selected=[args.unit]
            selection_requested=any(value is not None for value in (args.topology,args.assignment_stop,args.producer,args.placement,args.colour)) or args.assignment_start != 0
            if selection_requested:
                from neutrinomass.database.planning import generation_plan
                if args.unit is not None or args.stage not in {"generation","decoded","local-reduction"}:
                    parser.error("Topology/assignment selection requires a generation stage and no --unit")
                selected=[]
                for i,item in enumerate(generation_plan(args.operator,args.assignment_size)):
                    if any(value is not None and item[key] != value for key,value in (
                        ("topology",args.topology),("producer",args.producer),("placement",args.placement),("colour",args.colour))):
                        continue
                    if item["assignment_stop"] <= args.assignment_start or (args.assignment_stop is not None and item["assignment_start"] >= args.assignment_stop):
                        continue
                    if item["assignment_start"] < args.assignment_start or (args.assignment_stop is not None and item["assignment_stop"] > args.assignment_stop):
                        parser.error("Assignment limits must align with the selected --assignment-size plan boundaries")
                    selected.append(i)
                if not selected:
                    parser.error("Selection has no generation units")
            for count,unit in enumerate(selected):
                report = census_operator(
                    args.operator,args.historical_path,args.output_dir,
                    hash_seed=os.environ.get("PYTHONHASHSEED", ""),
                    stage=args.stage,unit=unit,assignment_size=args.assignment_size,
                    stage_unit=args.stage_unit,
                    resume=args.resume or count>0,budget=budget,scratch=args.scratch,allow_dirty=args.allow_dirty,
                )
    except BudgetExpired as error:
        print(json.dumps(dict(event="continuation_saved",reason=str(error))),flush=True)
        sys.exit(75)
    if args.report:
        write_report(args.report, report)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
