#!/usr/bin/env python3
"""Persistent census queue workers. Every command requires a Slurm allocation."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import time
import traceback

from neutrinomass.database.checkpoints import WorkBudget,BudgetExpired,atomic_json,digest_file
from neutrinomass.database.provenance import runtime_provenance,scientific_contract
from neutrinomass.database.planning import generation_plan
from neutrinomass.database.pipeline import run_operator
from neutrinomass.database.work_queue import WorkQueue,progress_signature


def prepare(config,provenance):
    queue=WorkQueue(config["queue"])
    historical=Path(config["historical"])
    for operator in config["operators"]:
        source=historical/f"op_{operator}.dat"
        if digest_file(source)!=config["historical_hashes"][operator]:
            raise ValueError(f"HISTORICAL_INPUT_HASH_MISMATCH: {operator}")
        plan=generation_plan(operator,config["assignment_size"])
        for seed in config["seeds"]:
            identifiers=[]
            for ordinal in range(len(plan)):
                identifier=f"{operator}/seed-{seed}/generation/{ordinal:08d}"
                queue.add(identifier,dict(operator=operator,unit=ordinal,stage="local-reduction"),
                    phase=config["phase"],seed=seed)
                identifiers.append(identifier)
            # Explicitly runnable finalisation; no ephemeral afterok chain.
            queue.add(f"{operator}/seed-{seed}/finalise",dict(operator=operator,unit=None,stage="all"),
                phase=config["phase"],seed=seed,dependencies=identifiers)
    atomic_json(Path(config["root"])/"prepared.json",dict(
        source_commit=provenance["source_commit"],config_sha256=config["_sha256"],
        operators=config["operators"],seeds=config["seeds"],queue=queue.summary()))


def work(config,provenance,*,seconds=5400):
    queue=WorkQueue(config["queue"])
    worker=os.environ["SLURM_JOB_ID"]
    seed=int(os.environ["PYTHONHASHSEED"])
    deadline=time.monotonic()+seconds
    contract=scientific_contract(provenance)
    scratch=Path(os.environ.get("SLURM_TMPDIR",os.environ.get("TMPDIR","/tmp")))/f"erratum-{worker}"
    scratch.mkdir(parents=True,exist_ok=True)
    while deadline-time.monotonic()>300:
        selected=queue.claim(worker,phase=config["phase"],seed=seed)
        if selected is None:
            break
        payload=selected["payload"]
        operator=payload["operator"]
        output=Path(config["root"])/f"seed-{seed}"/operator
        budget=WorkBudget(min(3300,deadline-time.monotonic()))
        budget.scratch_root=scratch
        try:
            with budget.signals():
                result=run_operator(operator,Path(config["historical"])/f"op_{operator}.dat",
                    output,contract,unit=payload["unit"],stage=payload["stage"],
                    assignment_size=config["assignment_size"],resume=(output/"plan.json").exists(),
                    budget=budget,scratch=scratch)
        except BudgetExpired:
            owned=(output if payload["unit"] is None else
                   output/"units"/f"{payload['unit']:08d}")
            status=queue.continuation(selected["id"],worker,progress_signature(owned))
            print(json.dumps(dict(event="saved_continuation",unit=selected["id"],status=status)),flush=True)
            return 75 if status=="pending" else 2
        except BaseException as error:
            queue.fail(selected["id"],worker,dict(code="WORKER_FAILURE",type=type(error).__name__,
                message=str(error),traceback=traceback.format_exc()))
            raise
        queue.finish(selected["id"],worker,dict(output=str(output),stage=payload["stage"],
            unit=payload["unit"],report_sha256=digest_file(output/"report.json")
            if payload["stage"]=="all" else None))
        print(json.dumps(dict(event="unit_complete",unit=selected["id"],
            elapsed_seconds=time.monotonic()-budget.started)),flush=True)
    return 0


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command",choices=("prepare","work"))
    parser.add_argument("config",type=Path)
    args=parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Spartan computations require Slurm")
    provenance=runtime_provenance()
    branch=subprocess.run(["git","symbolic-ref","-q","HEAD"],capture_output=True,text=True)
    if branch.returncode==0:
        parser.error("Worker checkout must be detached")
    config=json.loads(args.config.read_text())
    config["_sha256"]=digest_file(args.config)
    if config["source_commit"] != provenance["source_commit"]:
        parser.error("Config and tested source commit differ")
    if config["phase"] not in {"validation","production"}:
        parser.error("Unknown execution phase")
    if os.environ.get("SLURM_JOB_PARTITION")!="sapphire" or os.environ.get("SLURM_JOB_ACCOUNT")!="punim0011":
        parser.error("Workers require sapphire and account punim0011")
    if int(os.environ.get("SLURM_CPUS_PER_TASK","1"))!=1:
        parser.error("Workers require one allocated CPU")
    atomic_json(Path(config["root"])/"environments"/(os.environ["SLURM_JOB_ID"]+".json"),provenance)
    if args.command=="prepare":
        prepare(config,provenance)
        return 0
    return work(config,provenance)


if __name__=="__main__":
    raise SystemExit(main())
