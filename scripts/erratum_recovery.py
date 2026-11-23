"""Fault replicas and full-result substitution checks for the 71p pilot.

Execute with PYTHONPATH pointing to the tested producer worktree. This script
records its own validator hash independently of that producer revision.
"""
import argparse
import json
import os
from pathlib import Path
import shutil

from neutrinomass.database.checkpoints import WorkBudget,BudgetExpired,atomic_json,digest_file,CheckpointStore
from neutrinomass.database.provenance import runtime_provenance,scientific_contract
from neutrinomass.database.pipeline import run_operator,artifact_path
from neutrinomass.database.coverage import historical_coverage
from neutrinomass.database.resumable_reduction import lines_from_paths,input_contract
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


def substitute(config,reference,replica,recovered_generation,unit,budget,contract):
    # Copy only verified generation units into a fresh directory. Every global
    # reduction, amplitude, historical comparison and publication is rerun.
    if not (replica/"plan.json").exists():
        if replica.exists() and any(replica.iterdir()):
            raise ValueError("Substitution directory is not empty")
        replica.mkdir(parents=True,exist_ok=True)
        for source in sorted((reference/"units").iterdir()):
            chosen=(recovered_generation/"units"/source.name if source.name==f"{unit:08d}" else source)
            shutil.copytree(chosen,replica/"units"/source.name,copy_function=os.link)
        shutil.copyfile(reference/"plan.json",replica/"plan.json")
    operator=config.get("recovery_operator","71p")
    report=run_operator(operator,Path(config["historical"])/f"op_{operator}.dat",replica,contract,
        assignment_size=config["assignment_size"],resume=True,budget=budget,
        scratch=os.environ.get("SLURM_TMPDIR"))
    old=json.loads((reference/"report.json").read_text())
    if stable_report_data(old)!=stable_report_data(report):
        raise ValueError("SUBSTITUTION_CHANGED_SCIENTIFIC_REPORT")
    counts={}
    for key,relation in (("exact_classes",are_equivalent_completions),
                         ("global_lagrangians",equivalent_global_lagrangians)):
        counts[key]=ordered_exact(old["artifacts"][key]["path"],report["artifacts"][key]["path"],relation)
    return dict(kind="substitution",scientific_report_agreement=True,ordered_exact_classes=counts,
        original_report=dict(path=str(reference/"report.json"),sha256=digest_file(reference/"report.json")),
        substituted_report=dict(path=str(replica/"report.json"),sha256=digest_file(replica/"report.json")))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind",choices=("generation","coverage","substitution"))
    parser.add_argument("config",type=Path);parser.add_argument("reference",type=Path)
    parser.add_argument("replica",type=Path);parser.add_argument("--unit",type=int,required=True)
    parser.add_argument("--fault-after",type=int);parser.add_argument("--recovered-generation",type=Path)
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
    budget=WorkBudget(3300,fault_after=args.fault_after)
    try:
        with budget.signals():
            if args.kind=="generation":
                result=generation_replica(config,args.reference,args.replica,args.unit,budget,contract)
            elif args.kind=="coverage":
                result=coverage_replica(args.reference,args.replica,args.unit,budget,stage_contract)
            else:
                if args.recovered_generation is None:
                    parser.error("Substitution needs a verified recovered generation replica")
                result=substitute(config,args.reference,args.replica,args.recovered_generation,args.unit,budget,contract)
    except BudgetExpired:
        checkpoints=[dict(path=str(p),sha256=digest_file(p),progress=json.loads(p.read_text())["progress"])
                     for p in sorted(args.replica.rglob("checkpoint.json"))]
        atomic_json(args.replica/"interruptions"/(os.environ["SLURM_JOB_ID"]+".json"),dict(kind=args.kind,unit=args.unit,
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
