#!/usr/bin/env python3
"""Bounded queue workers for census verification, filtering and finalisation."""
import argparse
import json
import os
from pathlib import Path
import time
import traceback

from neutrinomass.database.checkpoints import CheckpointStore,WorkBudget,BudgetExpired,atomic_json,digest_file
from neutrinomass.database.campaign_finalisation import (
    completed_stage,SpeciesIndex,filter_species_unit,load_reports,
)
from neutrinomass.database.pipeline import artifact_path
from neutrinomass.database.provenance import runtime_provenance,scientific_contract
from neutrinomass.database.resumable_reduction import lines_from_paths
from neutrinomass.database.work_queue import WorkQueue,progress_signature
from scripts.finalise_erratum import run


def prepare(config):
    queue=WorkQueue(config["queue"])
    dependencies={"verify":(),"published":("verify",),"species-index":("verify",),
        "global-all":("verify",),"global-published":("published",),
        "comparison":("published",),"historical-replay":("published",)}
    for stage,parents in dependencies.items():
        queue.add(stage,dict(stage=stage),phase="finalisation",seed=0,dependencies=parents)
    atomic_json(Path(config["root"])/"prepared.json",dict(config_sha256=config["_sha256"],queue=queue.summary()))


def expand_species_queue(config,index_state,contract,*,budget):
    root=Path(config["root"]);queue=WorkQueue(config["queue"])
    work=artifact_path(root/"species-index",index_state,"work")
    identity=dict(contract,algorithm="species-queue-expansion-v1",work_sha256=digest_file(work))
    with CheckpointStore(root/"species-queue-expansion",identity) as store:
        if store.state["complete"]:return store.state
        progress=dict(store.state["progress"] or dict(added=0))
        identifiers=[]
        for ordinal,line in enumerate(lines_from_paths([work]),1):
            unit=json.loads(line);identifier="species/"+unit["identifier"]
            identifiers.append(identifier)
            if ordinal<=progress["added"]:continue
            queue.add(identifier,dict(stage="species-filter",work=unit),phase="finalisation",seed=0,
                      dependencies=("species-index",))
            progress["added"]=ordinal
            budget.boundary(lambda:store.save(progress),progress=dict(stage="species-queue-expansion",**progress))
        queue.add("species-collect",dict(stage="species-collect"),phase="finalisation",seed=0,
                  dependencies=("species-index",*identifiers))
        queue.add("global-species-aware",dict(stage="global-species-aware"),phase="finalisation",seed=0,
                  dependencies=("species-collect",))
        queue.add("report",dict(stage="report"),phase="finalisation",seed=0,
                  dependencies=("global-all","global-published","global-species-aware","comparison","historical-replay"))
        store.save(progress,complete=True)
        return store.state


def work(config,provenance,*,seconds=5400):
    queue=WorkQueue(config["queue"]);root=Path(config["root"])
    worker=os.environ["SLURM_JOB_ID"]
    contract=scientific_contract(provenance)
    manifest_path=Path(config["manifest"])
    if digest_file(manifest_path)!=config["manifest_sha256"]:
        raise ValueError("FINALISATION_MANIFEST_HASH_MISMATCH")
    manifest=json.loads(manifest_path.read_text())
    scratch=Path(os.environ.get("SLURM_TMPDIR",os.environ.get("TMPDIR","/tmp")))/f"erratum-finalise-{worker}"
    scratch.mkdir(parents=True,exist_ok=True)
    deadline=time.monotonic()+seconds
    index=None;scales=None
    try:
        while deadline-time.monotonic()>300:
            selected=queue.claim(worker,phase="finalisation",seed=0)
            if selected is None:break
            stage=selected["payload"]["stage"]
            budget=WorkBudget(min(3300,deadline-time.monotonic()));budget.scratch_root=scratch
            owned=root
            try:
                with budget.signals():
                    if stage=="species-filter":
                        if index is None:
                            state=completed_stage(root/"species-index",contract)
                            index=SpeciesIndex(dict(path=str(artifact_path(root/"species-index",state,"database")),
                                sha256=state["artifacts"]["database"]["sha256"])).__enter__()
                            reports=load_reports(manifest["entries"],verify_artifacts=False)
                            scales={op:seeds["0"]["operator_scale_tev"] for op,seeds in reports.items()}
                        unit=selected["payload"]["work"]
                        owned=root/"species-units"/unit["identifier"]
                        filter_species_unit(unit,index,scales,owned,contract,budget=budget)
                    else:
                        run(manifest,root,contract,stage=stage,budget=budget,scratch=scratch)
                        if stage=="species-index":
                            expand_species_queue(config,completed_stage(root/"species-index",contract),contract,budget=budget)
            except BudgetExpired:
                # Other tasks can advance concurrently. Only this task's own
                # checkpoints may demonstrate progress on its continuation.
                if stage!="species-filter":
                    locations={"verify":"seed-agreement","species-index":"species-index",
                        "species-collect":"species-results","global-all":"globals/all",
                        "global-published":"globals/published","global-species-aware":"globals/species-aware",
                        "comparison":"comparison","historical-replay":"historical-replay","report":"publication",
                        "published":"published"}
                    progress={locations[stage]:progress_signature(root/locations[stage])}
                    if stage=="species-index":progress["queue_expansion"]=progress_signature(root/"species-queue-expansion")
                else:progress=progress_signature(owned)
                status=queue.continuation(selected["id"],worker,progress)
                print(json.dumps(dict(event="saved_continuation",unit=selected["id"],status=status)),flush=True)
                return 75 if status=="pending" else 2
            except BaseException as error:
                queue.fail(selected["id"],worker,dict(code="FINALISATION_WORKER_FAILURE",
                    type=type(error).__name__,message=str(error),traceback=traceback.format_exc()))
                raise
            queue.finish(selected["id"],worker,dict(stage=stage,output=str(root)))
            print(json.dumps(dict(event="unit_complete",unit=selected["id"],elapsed_seconds=time.monotonic()-budget.started)),flush=True)
    finally:
        if index is not None:index.__exit__(None,None,None)
    return 0


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command",choices=("prepare","work"));parser.add_argument("config",type=Path)
    args=parser.parse_args()
    if (not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_JOB_PARTITION")!="sapphire"
            or os.environ.get("SLURM_JOB_ACCOUNT")!="punim0011"
            or int(os.environ.get("SLURM_CPUS_PER_TASK","1"))!=1
            or os.environ.get("PYTHONHASHSEED")!="0"):
        parser.error("Finalisation requires one sapphire CPU under punim0011, hash seed zero")
    provenance=runtime_provenance()
    config=json.loads(args.config.read_text());config["_sha256"]=digest_file(args.config)
    if config["source_commit"]!=provenance["source_commit"] or config["phase"]!="finalisation":
        parser.error("Finalisation config and pinned source differ")
    provenance["worker_driver_sha256"]=digest_file(__file__)
    atomic_json(Path(config["root"])/"environments"/(os.environ["SLURM_JOB_ID"]+".json"),provenance)
    if args.command=="prepare":prepare(config);return 0
    return work(config,provenance)


if __name__=="__main__":raise SystemExit(main())
