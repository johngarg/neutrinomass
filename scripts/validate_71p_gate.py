"""Verify the complete pilot before admitting any production allocation.

Run in Slurm with the tested producer on PYTHONPATH. Exact cross-seed checks
are themselves checkpointed, so this gate can be continued independently.
"""
import argparse
from collections import Counter
import json
import os
from pathlib import Path

from neutrinomass.database.campaign import CampaignLedger
from neutrinomass.database.checkpoints import CheckpointStore,WorkBudget,BudgetExpired,atomic_json,digest_file
from neutrinomass.database.coverage import historical_coverage
from neutrinomass.database.pipeline import artifact_path
from neutrinomass.database.provenance import runtime_provenance,scientific_contract
from neutrinomass.database.resumable_reduction import reduce_resumable
from neutrinomass.database.resumable_reduction import lines_from_paths
from neutrinomass.database.filtering import published_model_key
from neutrinomass.database.work_queue import WorkQueue
from rebuild_completion_database import stable_report_data


class GateFailure(ValueError):
    pass


def require(condition,code):
    if not condition:
        raise GateFailure(code)


def verify_reference(reference):
    path=Path(reference["path"])
    require(path.is_file() and digest_file(path)==reference["sha256"],"GATE_ARTIFACT_HASH_MISMATCH")
    return path


def report_gate(report,*,seed,producer,unit_count):
    require(report["operator"]=="71p" and str(report["hash_seed"])==str(seed),"GATE_OPERATOR_OR_SEED_MISMATCH")
    require(report["source_commit"]==producer,"GATE_PRODUCER_MISMATCH")
    require(report["stage_manifest"]["generation_units"]==unit_count,"GATE_GENERATION_UNITS_INCOMPLETE")
    require(len(report["stage_manifest"]["topologies"])==264,"GATE_TOPOLOGY_COVERAGE_INCOMPLETE")
    require(report["records"]["generator"]==758208,"GATE_RAW_REFERENCE_DEVIATION_REQUIRES_EXPLANATION")
    require(report["historical"]["missing"]==0 and not report["historical"]["amplitude_unsupported"],
            "GATE_HISTORICAL_OR_UNSUPPORTED_PHYSICS")
    require(report["historical"]["valid_classes"]==report["historical"]["reproduced"],"GATE_HISTORICAL_ACCOUNTING")
    require(report["exact_classes"]["all"]>0,"GATE_EMPTY_PHYSICAL_RESULT")


def verify_tree(root,producer):
    peak=0;stages=0
    for path in sorted(root.rglob("checkpoint.json")):
        data=json.loads(path.read_text())
        require(data["complete"],"GATE_INCOMPLETE_CHECKPOINT")
        require(data["contract"].get("source_commit")==producer,"GATE_CHECKPOINT_PRODUCER_MISMATCH")
        with CheckpointStore(path.parent,data["contract"]) as checked:
            peak=max(peak,checked.state.get("memory",{}).get("peak_rss_mib",0))
        stages+=1
    return dict(verified_stages=stages,peak_checkpoint_rss_mib=peak)


def substitution_compatibility_gate(data,producer,*,operator="71p"):
    """Validate the actual dependency bytes behind a downstream reuse claim."""
    proof_path=verify_reference(data["dependency_compatibility"])
    proof=json.loads(proof_path.read_text())
    identity=proof["contract"]
    require(proof["complete"] and identity.get("source_commit")==producer
            and identity.get("algorithm")=="recovery-substitution-setup-v2",
            "GATE_SUBSTITUTION_COMPATIBILITY_INCOMPLETE")
    with CheckpointStore(proof_path.parent,identity):pass
    require(identity["all_downstream_inputs_identical"],"GATE_SUBSTITUTION_INPUTS_CHANGED")
    require(identity["original_report_sha256"]==data["original_report"]["sha256"],
            "GATE_SUBSTITUTION_REFERENCE_CHANGED")
    original_root=Path(data["original_report"]["path"]).parent
    require(digest_file(original_root/"plan.json")==identity["original_plan_sha256"],
            "GATE_SUBSTITUTION_PLAN_CHANGED")
    generation={key:json.loads(verify_reference(ref).read_text())
                for key,ref in identity["generation_inputs"].items()}
    a,b=generation["original"],generation["recovered"]
    require(a["complete"] and b["complete"] and a["contract"]==b["contract"]
            and a["contract"]["source_commit"]==producer and a["contract"]["operator"]==operator
            and a["progress"]["records"]==b["progress"]["records"],
            "GATE_SUBSTITUTION_GENERATION_CHANGED")
    roots={key:Path(ref["path"]).parent.parent for key,ref in identity["generation_inputs"].items()}
    require(roots["original"].resolve()==(original_root/"units"/f'{identity["unit"]:08d}').resolve(),
            "GATE_SUBSTITUTION_UNIT_CHANGED")
    pairs=identity["consumer_inputs"]
    require(len(pairs)==3 and {(p["stage"],p["artifact"]) for p in pairs}==
            {("decoded","output"),("decoded","evidence"),("structural","output")},
            "GATE_SUBSTITUTION_CONSUMER_INPUTS_MISSING")
    for pair in pairs:
        paths={}
        for key in ("original","recovered"):
            path=verify_reference(pair[key]);stage=roots[key]/pair["stage"]
            manifest=json.loads((stage/"checkpoint.json").read_text())
            with CheckpointStore(stage,manifest["contract"]) as store:
                require(store.state["complete"] and artifact_path(stage,store.state,pair["artifact"]).resolve()==path.resolve(),
                        "GATE_SUBSTITUTION_CONSUMER_PATH_CHANGED")
            paths[key]=path
        with paths["original"].open("rb") as left,paths["recovered"].open("rb") as right:
            while True:
                x,y=left.read(1024**2),right.read(1024**2)
                require(x==y,"GATE_SUBSTITUTION_CONSUMER_BYTES_CHANGED")
                if not x:break
    for reference in identity["source_checkpoints"]:
        path=verify_reference(reference);manifest=json.loads(path.read_text())
        require(manifest["complete"] and manifest["contract"]["source_commit"]==producer,
                "GATE_SUBSTITUTION_SOURCE_STAGE_INCOMPLETE")
        with CheckpointStore(path.parent,manifest["contract"]):pass


def recovery_gate(directories,producer):
    found={};references=[]
    for directory in directories:
        path=Path(directory)/"recovery-evidence.json"
        data=json.loads(path.read_text())
        require(data["complete"] and data["producer"]["source_commit"]==producer,"GATE_RECOVERY_NOT_VERIFIED")
        require(data.get("job_id"),"GATE_RECOVERY_NOT_IN_SLURM")
        kind=data["kind"]
        if kind in {"generation","reduction","coverage"}:
            interruptions=[json.loads(p.read_text()) for p in (Path(directory)/"interruptions").glob("*.json")]
            injected=[item for item in interruptions if item.get("injected") and
                      item.get("method")=="SIGKILL" and item.get("child_returncode")==-9
                      and item.get("kind")==kind and item.get("unit")==data["unit"] and item.get("checkpoints")]
            require(injected,"GATE_NO_INJECTED_RECOVERY")
            for interruption in injected:
                for reference in interruption["checkpoints"]:
                    checkpoint=verify_reference(reference)
                    original=json.loads(checkpoint.read_text())
                    require(not original["complete"] and original["progress"]==reference["progress"],
                            "GATE_FAULT_CHECKPOINT_NOT_PARTIAL")
                    require(original["contract"]["source_commit"]==producer,"GATE_FAULT_PRODUCER_MISMATCH")
                    require(original["contract"].get("operator")=="71p","GATE_FAULT_OPERATOR_MISMATCH")
                    if kind=="reduction":
                        pending=original["progress"].get("pending")
                        require(pending and pending["cursor"][1]>=0
                                and pending["ordinal"]==original["progress"]["processed"]+1,
                                "GATE_REDUCTION_FAULT_NOT_INSIDE_BUCKET")
                    with CheckpointStore(checkpoint.parent,original["contract"]):pass
        if kind=="generation":
            require(data["ordered_exact_agreement"],"GATE_GENERATION_RECOVERY_DIFFERS")
        elif kind=="coverage":
            require(data["exact_output_agreement"] and data["full_coverage_substitution_agreement"],"GATE_COVERAGE_RECOVERY_DIFFERS")
        elif kind=="reduction":
            require(data["ordered_exact_agreement"] and data["records"]==data["progress"]["exact_classes"],
                    "GATE_REDUCTION_RECOVERY_DIFFERS")
            for key in ("original_output","recovered_output"):verify_reference(data[key])
        elif kind=="substitution":
            require(data["scientific_report_agreement"] and data["ordered_exact_classes"],"GATE_SUBSTITUTED_RESULT_DIFFERS")
            for key in ("original_report","substituted_report"):
                report=json.loads(verify_reference(data[key]).read_text())
                require(report["operator"]=="71p" and report["source_commit"]==producer,
                        "GATE_SUBSTITUTION_OPERATOR_MISMATCH")
                require(data["ordered_exact_classes"]["exact_classes"]==report["exact_classes"]["all"]
                        and data["ordered_exact_classes"]["global_lagrangians"]==report["global_lagrangians"],
                        "GATE_SUBSTITUTION_CLASS_COUNTS_DIFFER")
            require(data.get("method") in {"byte-identical-input-reuse","full-downstream-recomputation"},
                    "GATE_SUBSTITUTION_METHOD_UNSPECIFIED")
            if data["method"]=="byte-identical-input-reuse":substitution_compatibility_gate(data,producer)
        found[kind]=data
        references.append(dict(path=str(path),sha256=digest_file(path)))
    require(set(found)=={"generation","reduction","coverage","substitution"},"GATE_RECOVERY_EVIDENCE_INCOMPLETE")
    return references


def resource_gate(root,ledger,config_path):
    durations=[];scratch=[];capacities=[];memory=[]
    for path in sorted((root/"logs").glob("*.log")):
        for line in path.read_text(errors="replace").splitlines():
            if not line.startswith("{"):
                continue
            try:data=json.loads(line)
            except ValueError:continue
            if data.get("event")=="unit_complete":
                durations.append(data["elapsed_seconds"])
            if data.get("event")=="scratch_capacity":
                capacities.append(data["free_bytes"])
            if "scratch_bytes" in data:scratch.append(data["scratch_bytes"])
            if "peak_rss_mib" in data:memory.append(data["peak_rss_mib"])
    require(durations and capacities and memory,"GATE_MISSING_RESOURCE_MEASUREMENTS")
    require(max(durations)<3600,"GATE_WORK_UNIT_NEEDS_SUBDIVISION")
    require(max(memory)<6*1024,"GATE_MEMORY_HEADROOM_INSUFFICIENT")
    require(scratch and max(scratch)<0.8*min(capacities),"GATE_SCRATCH_HEADROOM_NOT_DEMONSTRATED")
    jobs=[j for j in ledger["jobs"].values()
          if j.get("queue_config")==str(Path(config_path).resolve()) and j.get("command") in {"prepare","work"}]
    # All allocations, including diagnostic failures, stay charged. A known
    # startup failure has a separate diagnosis and no scientific output.
    require(all("allocated_core_seconds" in j for j in jobs),"GATE_UNRECONCILED_ALLOCATIONS")
    require(all(j["allocated_core_seconds"]<=6000*j["cpus"] for j in jobs),"GATE_WALLTIME_HEADROOM_INSUFFICIENT")
    return dict(max_completed_unit_seconds=max(durations),peak_sampled_rss_mib=max(memory),
        peak_sampled_scratch_bytes=max(scratch),minimum_starting_scratch_free_bytes=min(capacities),
        completed_units_measured=len(durations))


def controls_gate(config):
    base=Path(config["root"]).parent
    producer=config["source_commit"]
    expected={"1":(3,3),"D3":(5,4),"D5b":(44,30),"D20":(279,129)}
    evidence=[];reports={}
    for seed in (0,1):
        complete=base/"validation"/producer/f"seed-{seed}-complete.json"
        data=json.loads(complete.read_text())
        require(data["suite_and_controls_complete"] and data["source_commit"]==producer,"GATE_FULL_SUITE_INCOMPLETE")
        environment=base/"validation"/producer/f"seed-{seed}-environment.json"
        env=json.loads(environment.read_text())
        require(env["source_commit"]==producer and not env["dirty_worktree"] and
                env["python_version"]=="3.11.3" and env["source_directory"]==config["source"],"GATE_VALIDATION_ENVIRONMENT_MISMATCH")
        log=base/"logs"/f"full-suite-{producer[:7]}-seed-{seed}-{env['slurm_job_id']}.log"
        require("307 passed" in log.read_text(),"GATE_FULL_SCIENTIFIC_SUITE_NOT_RECORDED")
        evidence.extend(dict(path=str(path),sha256=digest_file(path)) for path in (complete,environment,log))
        reports[seed]={}
        for operator,(physical,historical) in expected.items():
            path=base/"controls"/producer/f"seed-{seed}"/operator/"report.json"
            report=json.loads(path.read_text())
            require(report["source_commit"]==producer and report["exact_classes"]["all"]==physical and
                    report["historical"]["reproduced"]==historical and report["historical"]["missing"]==0,
                    "GATE_SMALL_OPERATOR_REGRESSION")
            for artifact in report["artifacts"].values():verify_reference(artifact)
            reports[seed][operator]=report
            evidence.append(dict(path=str(path),sha256=digest_file(path)))
    require(all(stable_report_data(reports[0][op])==stable_report_data(reports[1][op]) for op in expected),
            "GATE_SMALL_OPERATOR_SEED_DISAGREEMENT")
    return evidence


def validate(config_path,recoveries,baseline,*,budget):
    config=json.loads(Path(config_path).read_text());root=Path(config["root"])
    provenance=runtime_provenance()
    producer=config["source_commit"]
    require(provenance["source_commit"]==producer,"GATE_WRONG_SCIENTIFIC_VALIDATOR")
    summary=WorkQueue(config["queue"]).summary()
    require(not summary["blocked"] and all(row["status"]=="complete" for row in summary["counts"]),"GATE_WORK_QUEUE_INCOMPLETE")
    historical=json.loads(Path(baseline).read_text())
    require(historical["passed"] and not historical["per_operator_differences"],"GATE_HISTORICAL_BASELINE_FAILED")
    require((historical["literal_unfiltered_models"],historical["filtered_models"],
             historical["historical_interaction_multisets"],historical["historical_filtered_first_representative_lagrangians"])
            ==(141989,11216,430810,11483),"GATE_PUBLISHED_BASELINE_COUNTS")
    controls=controls_gate(config)
    reports=[];verified=[]
    directory=root/"acceptance"
    for seed in (0,1):
        operator_root=root/f"seed-{seed}"/"71p"
        report=json.loads((operator_root/"report.json").read_text())
        plan=json.loads((operator_root/"plan.json").read_text())
        report_gate(report,seed=seed,producer=producer,unit_count=len(plan["units"]))
        require(set(u["topology"] for u in plan["units"])==set(report["stage_manifest"]["topologies"]),"GATE_PLAN_TOPOLOGIES_DIFFER")
        for ordinal,selection in enumerate(plan["units"]):
            for stage in ("generation","decoded","structural"):
                unit=operator_root/"units"/f"{ordinal:08d}"/stage/"checkpoint.json"
                require(unit.is_file(),"GATE_REQUIRED_UNIT_CHECKPOINT_MISSING")
                data=json.loads(unit.read_text())
                require(data["complete"],"GATE_REQUIRED_UNIT_NOT_COMPLETE")
                if stage=="generation":
                    require(data["contract"]["selection"]==selection,"GATE_UNIT_ASSIGNMENT_RANGE_CHANGED")
        for artifact in report["artifacts"].values():verify_reference(artifact)
        verified.append(verify_tree(operator_root,producer));reports.append(report)
    require(stable_report_data(reports[0])==stable_report_data(reports[1]),"GATE_SEED_REPORT_DISAGREEMENT")
    model_sets=[{published_model_key(json.loads(line)["fields"])
                 for line in lines_from_paths([r["artifacts"]["democratic_models"]["path"]])} for r in reports]
    require(model_sets[0]==model_sets[1],"GATE_SEED_MODEL_SETS_DIFFER")
    contract=dict(scientific_contract(provenance),validation="71p-exact-cross-seed-v1")
    exact=[]
    for source,target in ((0,1),(1,0)):
        coverage=historical_coverage([reports[target]["artifacts"]["exact_classes"]["path"]],
            [reports[source]["artifacts"]["exact_classes"]["path"]],directory/f"seed-{source}-to-{target}",
            contract,budget=budget,scratch=os.environ.get("SLURM_TMPDIR"))
        require(coverage["progress"]["counts"]=={"reproduced":reports[source]["exact_classes"]["all"]},"GATE_EXACT_SEED_CLASSES_DIFFER")
        exact.append(coverage["progress"])
    global_union=reduce_resumable([r["artifacts"]["global_lagrangians"]["path"] for r in reports],
        directory/"global-seed-union",contract,relation="global-lagrangian",budget=budget,scratch=os.environ.get("SLURM_TMPDIR"))
    require(global_union["progress"]["exact_classes"]==reports[0]["global_lagrangians"]==reports[1]["global_lagrangians"],"GATE_GLOBAL_SEED_CLASSES_DIFFER")
    recovery=recovery_gate(recoveries,producer)
    ledger=CampaignLedger(config["ledger"])
    accounting=ledger.summary()
    require(accounting["validation_charged_core_hours"]<=200,"GATE_VALIDATION_BUDGET_EXHAUSTED")
    resources=resource_gate(root,accounting,config_path)
    require(max(x["peak_checkpoint_rss_mib"] for x in verified)<6*1024,"GATE_CHECKPOINT_MEMORY_HEADROOM")
    result=dict(schema_version=1,passed=True,producer_commit=producer,validator_sha256=digest_file(__file__),
        verified=verified,exact_seed_coverage=exact,global_seed_union=global_union["progress"],resources=resources,
        recovery=recovery,controls=controls,baseline=dict(path=str(baseline),sha256=digest_file(baseline)),
        reports=[dict(path=str(root/f"seed-{seed}"/"71p/report.json"),sha256=digest_file(root/f"seed-{seed}"/"71p/report.json")) for seed in (0,1)],
        allocated_validation_core_hours=accounting["validation_charged_core_hours"],low_scale_genuineness_proved=False)
    path=directory/"gate.json";atomic_json(path,result)
    with ledger.transaction() as state:
        state["production_gate"]=dict(path=str(path),sha256=digest_file(path),producer_commit=producer)
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config",type=Path);parser.add_argument("--recovery",type=Path,action="append",required=True)
    parser.add_argument("--baseline",type=Path,required=True)
    args=parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):parser.error("The exact gate checks require Slurm")
    budget=WorkBudget(3300)
    try:
        with budget.signals():result=validate(args.config,args.recovery,args.baseline,budget=budget)
    except BudgetExpired:
        print(json.dumps(dict(event="gate_continuation_saved")),flush=True)
        raise SystemExit(75)
    print(json.dumps(result,sort_keys=True,indent=2))
