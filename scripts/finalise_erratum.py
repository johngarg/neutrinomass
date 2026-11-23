"""Resume verified census finalisation independently of any Slurm afterok chain."""
import argparse
import json
import os
from pathlib import Path

from neutrinomass.database.campaign_finalisation import (
    load_reports,verify_seeds,published_filter,build_species_index,SpeciesIndex,
    filter_species_unit,collect_species_results,global_reduction,model_comparison,
    historical_counting_replay,completed_stage,publish_summary,verify_reference,
)
from neutrinomass.database.checkpoints import (
    CheckpointStore,WorkBudget,BudgetExpired,atomic_json,digest_file,canonical_json,
)
from neutrinomass.database.pipeline import artifact_path
from neutrinomass.database.provenance import runtime_provenance,scientific_contract
from neutrinomass.database.resumable_reduction import lines_from_paths

STAGES=("verify","published","species-index","species-filter","species-collect",
        "global-all","global-published","global-species-aware","comparison","historical-replay","report")


def run(manifest,directory,contract,*,stage="all",unit=None,budget=None,scratch=None,expected=None):
    directory=Path(directory);budget=budget or WorkBudget()
    if unit is not None and stage!="species-filter":
        raise ValueError("FILTER_UNIT_SELECTION_REQUIRES_SPECIES_FILTER_STAGE")
    identity=dict(contract,algorithm="campaign-input-verification-v1",entries=manifest["entries"])
    verification=directory/"input-verification"
    if stage in {"verify","all"}:
        reports=load_reports(manifest["entries"],expected=expected)
        with CheckpointStore(verification,identity) as store:
            if not store.state["complete"]:
                proof=verify_seeds(reports,directory/"seed-agreement",contract,budget=budget,scratch=scratch)
                output=store.publish_lines(canonical_json(row) for row in proof)
                store.save(dict(operators=len(reports),seeds=[0,1]),artifacts={"output":output},complete=True)
    else:
        with CheckpointStore(verification,identity) as store:
            if not store.state["complete"]:raise ValueError("CAMPAIGN_INPUT_VERIFICATION_REQUIRED")
        reports=load_reports(manifest["entries"],expected=expected,verify_artifacts=False)
    if any(report["source_commit"]!=manifest["producer_commit"]
           for seeds in reports.values() for report in seeds.values()):
        raise ValueError("CAMPAIGN_INPUT_PRODUCER_MISMATCH")
    if stage=="verify":return dict(complete=True,stage=stage,operators=len(reports))
    scales={op:seeds["0"]["operator_scale_tev"] for op,seeds in reports.items()}
    if stage in {"published","all"}:
        published_filter(reports,directory/"published",contract)
    if stage in {"species-index","all"}:
        build_species_index(reports,directory/"species-index",contract,budget=budget,scratch=scratch)
    if stage in {"species-filter","all"}:
        index=completed_stage(directory/"species-index",contract)
        reference=dict(path=str(artifact_path(directory/"species-index",index,"database")),
                       sha256=index["artifacts"]["database"]["sha256"])
        selected=0
        with SpeciesIndex(reference) as verified:
            for ordinal,line in enumerate(lines_from_paths([artifact_path(directory/"species-index",index,"work")])):
                if unit is not None and ordinal!=unit:continue
                work=json.loads(line)
                filter_species_unit(work,verified,scales,directory/"species-units"/work["identifier"],contract,budget=budget)
                selected+=1
                # The unit's immutable output is already durable.
                budget.boundary(lambda:None,progress=dict(stage="species-filter",unit=ordinal))
        if unit is not None and selected!=1:raise ValueError("UNKNOWN_SPECIES_FILTER_UNIT")
    if stage in {"species-collect","all"}:
        index=completed_stage(directory/"species-index",contract)
        collect_species_results(directory/"species-index",index,directory/"species-units",directory/"species-results",
                                contract,scales=scales,budget=budget,scratch=scratch)
    for mode in ("all","published","species-aware"):
        if stage not in {"global-"+mode,"all"}:continue
        filter_path=None
        if mode!="all":
            root=directory/("published" if mode=="published" else "species-results")
            filter_path=artifact_path(root,completed_stage(root,contract))
        global_reduction(reports,directory/"globals"/mode,contract,mode=mode,
                         filter_path=filter_path,budget=budget,scratch=scratch)
    if stage in {"comparison","historical-replay","all"}:
        published=artifact_path(directory/"published",completed_stage(directory/"published",contract))
        if stage in {"comparison","all"}:
            model_comparison(reports,published,verify_reference(manifest["baseline"]),
                verify_reference(manifest["baseline_models"]),directory/"comparison",contract,budget=budget,scratch=scratch)
        if stage in {"historical-replay","all"}:
            historical_counting_replay(reports,published,directory/"historical-replay",contract,budget=budget,scratch=scratch)
    if stage in {"report","all"}:
        return publish_summary(manifest["entries"],directory,contract,manifest["baseline"],
                               expected=expected,budget=budget,scratch=scratch)
    return dict(complete=True,stage=stage,unit=unit)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest",type=Path);parser.add_argument("output",type=Path)
    parser.add_argument("--stage",choices=("all",*STAGES),default="all")
    parser.add_argument("--unit",type=int);parser.add_argument("--time-budget",type=float,default=3300)
    args=parser.parse_args()
    if (not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_JOB_PARTITION")!="sapphire"
            or os.environ.get("SLURM_JOB_ACCOUNT")!="punim0011"):
        parser.error("Census finalisation requires a sapphire Slurm allocation under punim0011")
    if not 0<args.time_budget<=5400:parser.error("Time budget must preserve the 90-minute exit boundary")
    provenance=runtime_provenance()
    provenance["finaliser_sha256"]=digest_file(__file__)
    atomic_json(args.output/"environments"/(os.environ["SLURM_JOB_ID"]+".json"),provenance)
    budget=WorkBudget(args.time_budget)
    scratch=os.environ.get("SLURM_TMPDIR")
    budget.scratch_root=scratch
    try:
        with budget.signals():
            result=run(json.loads(args.manifest.read_text()),args.output,scientific_contract(provenance),
                       stage=args.stage,unit=args.unit,budget=budget,scratch=scratch)
    except BudgetExpired:return 75
    print(json.dumps({key:result[key] for key in ("complete","complete_census","stage","unit","operator_count")
                      if key in result},sort_keys=True))
    return 0


if __name__=="__main__":raise SystemExit(main())
