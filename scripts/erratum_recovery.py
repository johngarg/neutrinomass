"""Fault replicas and full-result substitution checks for the 71p pilot.

Execute with PYTHONPATH pointing to the tested producer worktree. This script
records its own validator hash independently of that producer revision.
"""
import argparse
import json
import os
from pathlib import Path
import time

from neutrinomass.database.checkpoints import WorkBudget,BudgetExpired,atomic_json,digest_file,CheckpointStore
from neutrinomass.database.provenance import runtime_provenance,scientific_contract
from neutrinomass.database.pipeline import run_operator,artifact_path
from neutrinomass.database.coverage import historical_coverage
from neutrinomass.database.resumable_reduction import lines_from_paths,input_contract,reduce_resumable
from neutrinomass.database.serialization import loads_completion
from neutrinomass.completions.completions import are_equivalent_completions
from neutrinomass.database.global_equivalence import equivalent_global_lagrangians
from rebuild_completion_database import stable_report_data


def state(directory):
    return json.loads((Path(directory)/"checkpoint.json").read_text())


def ordered_exact(left,right,relation=are_equivalent_completions):
    from itertools import zip_longest
    count=0
    for a,b in zip_longest(lines_from_paths([left]),lines_from_paths([right])):
        if a is None or b is None or not relation(loads_completion(a),loads_completion(b)):
            raise ValueError("RECOVERED_ORDERED_CLASSES_DIFFER")
        count+=1
    return count


def bytes_equal(left,right):
    """Actual byte comparison, not a scientific fingerprint equality test."""
    with Path(left).open("rb") as a,Path(right).open("rb") as b:
        while True:
            x,y=a.read(1024**2),b.read(1024**2)
            if x!=y:return False
            if not x:return True


def substitution_setup(reference,replica,recovered,unit,budget,contract):
    selected=f"{unit:08d}"
    old=reference/"units"/selected;new=recovered/"units"/selected
    pairs=[]
    for stage,keys in (("decoded",("output","evidence")),("structural",("output",))):
        a,b=state(old/stage),state(new/stage)
        for key in keys:
            left=artifact_path(old/stage,a,key);right=artifact_path(new/stage,b,key)
            pairs.append(dict(stage=stage,artifact=key,
                original=dict(path=str(left),sha256=digest_file(left)),
                recovered=dict(path=str(right),sha256=digest_file(right)),
                bytes_identical=bytes_equal(left,right)))
    a,b=state(old/"generation"),state(new/"generation")
    generation_compatible=(a["complete"] and b["complete"] and a["contract"]==b["contract"]
                           and a["progress"]["records"]==b["progress"]["records"])
    identical=generation_compatible and all(pair["bytes_identical"] for pair in pairs)
    # When the actual consumer inputs are byte-identical, retain the original
    # downstream producer and validate every immutable artifact on reuse. This
    # avoids paying for an unnecessary third complete late-stage calculation.
    # If this stronger condition fails, recompute every downstream stage.
    sources=[]
    for checkpoint in sorted(reference.rglob("checkpoint.json")):
        relative=checkpoint.parent.relative_to(reference)
        if not identical and relative.parts[0]!="units":continue
        source=(recovered/relative if relative.parts[:2]==("units",selected) else checkpoint.parent)
        sources.append((relative,source))
    identity=dict(contract,algorithm="recovery-substitution-setup-v2",unit=unit,
        original_report_sha256=digest_file(reference/"report.json"),
        original_plan_sha256=digest_file(reference/"plan.json"),
        all_downstream_inputs_identical=identical,consumer_inputs=pairs,
        generation_inputs={key:dict(path=str(path/"checkpoint.json"),sha256=digest_file(path/"checkpoint.json"))
            for key,path in (("original",old/"generation"),("recovered",new/"generation"))},
        source_checkpoints=[dict(relative=str(relative),path=str(source/"checkpoint.json"),
            sha256=digest_file(source/"checkpoint.json")) for relative,source in sources])
    with CheckpointStore(replica/"substitution-setup",identity) as setup:
        if not setup.state["complete"]:
            progress=dict(setup.state["progress"] or dict(copied_stages=0))
            # Certify the identity before any copied output can survive a kill.
            # Repeating an unrecorded copy is safe: its bytes are checked below.
            if not (setup.directory/"checkpoint.json").exists():setup.save(progress)
            for ordinal,(relative,source) in enumerate(sources,1):
                if ordinal<=progress["copied_stages"]:continue
                original=state(source)
                with CheckpointStore(source,original["contract"]) as verified:
                    if not verified.state["complete"]:raise ValueError("SUBSTITUTION_SOURCE_STAGE_INCOMPLETE")
                    target=replica/relative;target.mkdir(parents=True,exist_ok=True)
                    for artifact in verified.state["artifacts"].values():
                        before=source/artifact["path"];after=target/artifact["path"]
                        if after.exists():
                            if not bytes_equal(before,after):raise ValueError("SUBSTITUTION_PARTIAL_COPY_CORRUPT")
                        else:os.link(before,after)
                    # Independent locks/manifests; immutable payloads can be
                    # hardlinked on GPFS without changing the original run.
                    atomic_json(target/"checkpoint.json",verified.state)
                progress["copied_stages"]=ordinal
                budget.boundary(lambda:setup.save(progress),progress=dict(stage="substitution-setup",**progress))
            atomic_json(replica/"plan.json",json.loads((reference/"plan.json").read_text()))
            setup.save(progress,complete=True)
    return identity


def generation_replica(config,reference,replica,unit,budget,contract):
    operator=config.get("recovery_operator","71p")
    run_operator(operator,Path(config["historical"])/f"op_{operator}.dat",replica,contract,
        assignment_size=config["assignment_size"],unit=unit,stage="local-reduction",
        resume=(replica/"plan.json").exists(),budget=budget,
        scratch=os.environ.get("SLURM_TMPDIR"))
    comparisons={}
    for stage in ("decoded","structural"):
        old=reference/"units"/f"{unit:08d}"/stage
        new=replica/"units"/f"{unit:08d}"/stage
        old_path=artifact_path(old,state(old));new_path=artifact_path(new,state(new))
        comparisons[stage]=dict(records=ordered_exact(old_path,new_path),
            original_sha256=digest_file(old_path),recovered_sha256=digest_file(new_path))
    return dict(kind="generation",unit=unit,ordered_exact_agreement=True,comparisons=comparisons)


def coverage_replica(reference,replica,unit,budget,contract):
    source=reference/"coverage"/"shards"/f"{unit:08d}"
    original=state(source)
    start,stop=original["contract"]["record_range"]
    split=state(reference/"coverage/inputs")
    chunk=artifact_path(reference/"coverage/inputs",split,f"chunk-{unit:08d}")
    report=json.loads((reference/"report.json").read_text())
    index=reference/"coverage/index"
    completed=historical_coverage([report["artifacts"]["exact_classes"]["path"]],[chunk],replica,contract,
        record_start=start,record_stop=stop,ordinal_offset=start,index_path=artifact_path(index,state(index),"database"),
        budget=budget,scratch=os.environ.get("SLURM_TMPDIR"))
    a=artifact_path(source,original);b=artifact_path(replica,completed)
    if list(lines_from_paths([a]))!=list(lines_from_paths([b])) or original["progress"]!=completed["progress"]:
        raise ValueError("RECOVERED_COVERAGE_DIFFERS")
    replacement=[]
    for directory in sorted((reference/"coverage/shards").iterdir()):
        replacement.append(b if directory.name==f"{unit:08d}" else artifact_path(directory,state(directory)))
    merge=replica/"substituted-coverage"
    with CheckpointStore(merge,dict(contract,algorithm="coverage-substitution-v1",inputs=input_contract(replacement))) as store:
        if not store.state["complete"]:
            output=store.publish_lines(lines_from_paths(replacement))
            store.save(dict(substituted_unit=unit),artifacts={"output":output},complete=True)
        substituted=artifact_path(merge,store.state)
    original_merged=report["artifacts"]["historical_coverage"]["path"]
    from itertools import zip_longest
    if any(a!=b for a,b in zip_longest(lines_from_paths([original_merged]),lines_from_paths([substituted]))):
        raise ValueError("SUBSTITUTED_COVERAGE_PUBLICATION_DIFFERS")
    return dict(kind="coverage",unit=unit,exact_output_agreement=True,
        original_sha256=digest_file(a),recovered_sha256=digest_file(b),progress=completed["progress"],
        full_coverage_substitution_agreement=True,substituted_coverage_sha256=digest_file(substituted))


def reduction_replica(reference,replica,budget,contract):
    """Repeat the actual global structural merge, including its bucket cursor."""
    plan=json.loads((reference/"plan.json").read_text())
    source=reference/"structural";original=state(source)
    paths=[]
    for unit in range(len(plan["units"])):
        stage=reference/"units"/f"{unit:08d}"/"structural"
        paths.append(artifact_path(stage,state(stage)))
    if original["contract"].get("ranking") is not None:
        raise ValueError("REDUCTION_REPLICA_REQUIRES_UNRANKED_PILOT")
    expected=dict(contract,algorithm="exact-reduction-v1",relation="physical",inputs=input_contract(paths),
                  ranking=None,maximum_rank=None)
    if original["contract"]!=expected:raise ValueError("REDUCTION_REPLICA_INPUT_CONTRACT_CHANGED")
    completed=reduce_resumable(paths,replica,contract,budget=budget,scratch=os.environ.get("SLURM_TMPDIR"))
    a=artifact_path(source,original);b=artifact_path(replica,completed)
    identical=bytes_equal(a,b)
    count=sum(1 for _ in lines_from_paths([b])) if identical else ordered_exact(a,b)
    if original["progress"]!=completed["progress"]:raise ValueError("RECOVERED_REDUCTION_PROGRESS_DIFFERS")
    return dict(kind="reduction",unit=0,ordered_exact_agreement=True,byte_identical=identical,
        records=count,progress=completed["progress"],
        original_output=dict(path=str(a),sha256=digest_file(a)),
        recovered_output=dict(path=str(b),sha256=digest_file(b)))


def substitute(config,reference,replica,recovered_generation,unit,budget,contract):
    if replica.resolve() in {reference.resolve(),recovered_generation.resolve()}:
        raise ValueError("SUBSTITUTION_REPLICA_OVERLAPS_SOURCE")
    if replica.exists() and any(replica.iterdir()) and not (replica/"substitution-setup/checkpoint.json").exists():
        # A kill between creating the setup lock and its initial manifest may
        # leave only this lock (or an incomplete atomic-publication temporary).
        files=list(replica.rglob("*"))
        safe=all(p==replica/"substitution-setup" or
                 (p.parent==replica/"substitution-setup" and p.is_file() and
                  (p.name=="worker.lock" or p.name.startswith(".publish-"))) for p in files)
        if not safe:raise ValueError("Substitution directory is not empty or has no resumable setup")
    compatibility=substitution_setup(reference,replica,recovered_generation,unit,budget,contract)
    operator=config.get("recovery_operator","71p")
    report=run_operator(operator,Path(config["historical"])/f"op_{operator}.dat",replica,contract,
        assignment_size=config["assignment_size"],resume=True,budget=budget,
        scratch=os.environ.get("SLURM_TMPDIR"))
    old=json.loads((reference/"report.json").read_text())
    if stable_report_data(old)!=stable_report_data(report):
        raise ValueError("SUBSTITUTION_CHANGED_SCIENTIFIC_REPORT")
    if compatibility["all_downstream_inputs_identical"]:
        # The original publication remains immutable. Give this validation
        # report explicit references to the substituted replica's actual bytes.
        report=json.loads(json.dumps(report))
        for artifact in (*report["artifacts"].values(),report["stage_manifest"]["plan"]):
            artifact["path"]=str(replica/Path(artifact["path"]).relative_to(reference))
        report["recovery_substitution"]=dict(original_report=str(reference/"report.json"),
            validator_sha256=digest_file(__file__),retained_producer=report["source_commit"])
    counts={}
    for key,relation in (("exact_classes",are_equivalent_completions),
                         ("global_lagrangians",equivalent_global_lagrangians)):
        original=Path(old["artifacts"][key]["path"])
        if compatibility["all_downstream_inputs_identical"]:
            copied=replica/original.relative_to(reference)
            if not bytes_equal(original,copied):raise ValueError("SUBSTITUTED_OUTPUT_BYTES_CHANGED")
            counts[key]=sum(1 for _ in lines_from_paths([copied]))
        else:
            counts[key]=ordered_exact(original,report["artifacts"][key]["path"],relation)
    publication=dict(contract,algorithm="substitution-validation-publication-v1",
        report=report,compatibility_sha256=digest_file(replica/"substitution-setup/checkpoint.json"))
    with CheckpointStore(replica/"substitution-publication",publication) as store:
        if not store.state["complete"]:
            output=store.publish_lines([json.dumps(report,sort_keys=True)])
            store.save(dict(ordered_exact_classes=counts),artifacts={"report":output},complete=True)
    atomic_json(replica/"report.json",report)
    return dict(kind="substitution",scientific_report_agreement=True,ordered_exact_classes=counts,
        method="byte-identical-input-reuse" if compatibility["all_downstream_inputs_identical"] else "full-downstream-recomputation",
        dependency_compatibility=dict(path=str(replica/"substitution-setup/checkpoint.json"),
            sha256=digest_file(replica/"substitution-setup/checkpoint.json")),
        original_report=dict(path=str(reference/"report.json"),sha256=digest_file(reference/"report.json")),
        substituted_report=dict(path=str(replica/"report.json"),sha256=digest_file(replica/"report.json")))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind",choices=("generation","reduction","coverage","substitution"))
    parser.add_argument("config",type=Path);parser.add_argument("reference",type=Path)
    parser.add_argument("replica",type=Path);parser.add_argument("--unit",type=int,required=True)
    parser.add_argument("--fault-after",type=int);parser.add_argument("--recovered-generation",type=Path)
    parser.add_argument("--checkpoint-seconds",type=float,default=300)
    parser.add_argument("--pause-after-pending-checkpoint",action="store_true",help=argparse.SUPPRESS)
    args=parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Recovery acceptance must run in Slurm")
    config=json.loads(args.config.read_text())
    provenance=runtime_provenance()
    if provenance["source_commit"]!=config["source_commit"]:
        parser.error("PYTHONPATH must select the tested producer revision")
    contract=scientific_contract(provenance)
    # Pipeline adds these operator-specific fields before calling its stages.
    operator=config.get("recovery_operator","71p")
    stage_contract=dict(contract,operator=operator,historical_sha256=digest_file(Path(config["historical"])/f"op_{operator}.dat"))
    budget=WorkBudget(3300,fault_after=args.fault_after,checkpoint_seconds=args.checkpoint_seconds)
    if args.pause_after_pending_checkpoint:
        if args.kind!="reduction":parser.error("Pending-checkpoint instrumentation only applies to reduction")
        snapshot=CheckpointStore.snapshot
        def observed_snapshot(store,database,progress,**kwargs):
            result=snapshot(store,database,progress,**kwargs)
            pending=progress.get("pending")
            if pending and pending["cursor"][1]>=0:
                # Fault instrumentation only: the durable checkpoint is already
                # published. Give the parent a deterministic SIGKILL window.
                time.sleep(2)
            return result
        CheckpointStore.snapshot=observed_snapshot
    try:
        with budget.signals():
            if args.kind=="generation":
                result=generation_replica(config,args.reference,args.replica,args.unit,budget,contract)
            elif args.kind=="coverage":
                result=coverage_replica(args.reference,args.replica,args.unit,budget,stage_contract)
            elif args.kind=="reduction":
                if args.unit!=0:parser.error("The global structural merge has unit index zero")
                result=reduction_replica(args.reference,args.replica,budget,stage_contract)
            else:
                if args.recovered_generation is None:
                    parser.error("Substitution needs a verified recovered generation replica")
                result=substitute(config,args.reference,args.replica,args.recovered_generation,args.unit,budget,contract)
    except BudgetExpired:
        checkpoints=[dict(path=str(p),sha256=digest_file(p),progress=json.loads(p.read_text())["progress"])
                     for p in sorted(args.replica.rglob("checkpoint.json"))]
        # Keep a parent's earlier SIGKILL certificate in the same allocation.
        atomic_json(args.replica/"interruptions"/("continuation-"+os.environ["SLURM_JOB_ID"]+".json"),dict(kind=args.kind,unit=args.unit,
            injected=args.fault_after is not None,job_id=os.environ["SLURM_JOB_ID"],checkpoints=checkpoints,
            validator_sha256=digest_file(__file__),producer=provenance))
        return 75
    result.update(validator_sha256=digest_file(__file__),producer=provenance,
                  job_id=os.environ["SLURM_JOB_ID"],complete=True)
    atomic_json(args.replica/"recovery-evidence.json",result)
    print(json.dumps(result,sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
