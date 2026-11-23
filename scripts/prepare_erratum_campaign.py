"""Create the gated execution/reuse manifest and seal completed census inputs.

This performs verification in Slurm. It never submits jobs; the bounded
controller consumes the resulting production configuration after the gate.
"""
import argparse
import json
import os
from pathlib import Path

from neutrinomass.database.checkpoints import CheckpointStore,atomic_json,digest_file
from neutrinomass.database.provenance import runtime_provenance,scientific_contract
from neutrinomass.database.rebuild import operator_registry
from neutrinomass.database.work_queue import WorkQueue


STAGES=("generation","decoded","local-reduction","structural-reduction","amplitudes",
        "historical-normalization","historical-reduction","historical-coverage",
        "global-lagrangians","packaging","publication")


def reference(path):
    path=Path(path).resolve()
    return dict(path=str(path),sha256=digest_file(path))


def verify(ref):
    path=Path(ref["path"])
    if not path.is_file() or digest_file(path)!=ref["sha256"]:
        raise ValueError("EXECUTION_MANIFEST_REFERENCE_CHANGED: "+str(path))
    return path


def gate_reference(config):
    ledger=json.loads(Path(config["ledger"]).read_text())
    ref=ledger.get("production_gate")
    if not ref:raise ValueError("EXECUTION_MANIFEST_REQUIRES_PASSED_71P_GATE")
    gate=json.loads(verify(ref).read_text())
    if not gate.get("passed") or gate["producer_commit"]!=config["source_commit"]:
        raise ValueError("EXECUTION_MANIFEST_GATE_PRODUCER_MISMATCH")
    return ref,gate


def verify_reuse(report_path,operator,seed,producer,historical_hash,validator,gate):
    ref=reference(report_path)
    approved=gate["reports"]+gate["controls"]
    if not any(item["path"]==ref["path"] and item["sha256"]==ref["sha256"] for item in approved):
        raise ValueError("REUSE_REPORT_NOT_VALIDATED_BY_PILOT_GATE")
    report=json.loads(Path(report_path).read_text())
    if (report["source_commit"]!=producer or report["operator"]!=operator
            or str(report["hash_seed"])!=str(seed) or report["schema_version"]!=2
            or report["scale_unit"]!="TeV" or report["historical"]["missing"]
            or report["historical"]["amplitude_unsupported"]):
        raise ValueError("REUSE_REPORT_SCIENTIFIC_IDENTITY_MISMATCH")
    if report["historical"]["artifact"]["sha256"]!=historical_hash:
        raise ValueError("REUSE_HISTORICAL_INPUT_CHANGED")
    verify(report["historical"]["artifact"])
    for artifact in report["artifacts"].values():verify(artifact)
    plan=json.loads(verify(report["stage_manifest"]["plan"]).read_text())
    base=plan["contract"]
    if base!=report["stage_manifest"]["contract"] or base["source_commit"]!=producer:
        raise ValueError("REUSE_STAGE_CONTRACT_CHANGED")
    root=Path(report_path).parent
    for unit in range(len(plan["units"])):
        for stage in ("generation","decoded","structural"):
            if not (root/"units"/f"{unit:08d}"/stage/"checkpoint.json").is_file():
                raise ValueError("REUSE_REQUIRED_UNIT_MISSING")
    stages=[]
    for checkpoint in sorted(root.rglob("checkpoint.json")):
        data=json.loads(checkpoint.read_text())
        if not data["complete"] or any(data["contract"].get(k)!=v for k,v in base.items()):
            raise ValueError("REUSE_STAGE_INCOMPATIBLE")
        with CheckpointStore(checkpoint.parent,data["contract"]):pass
        stages.append(dict(stage=str(checkpoint.parent.relative_to(root)),checkpoint=reference(checkpoint),
            producer=base,validator=validator,reuse=True,code="IDENTICAL_PRODUCER_AND_VERIFIED_FULL_RUN"))
    if not stages:raise ValueError("REUSE_REQUIRES_UNDERLYING_STAGES")
    if len(plan["units"])!=report["stage_manifest"]["generation_units"]:
        raise ValueError("REUSE_PLAN_UNIT_COUNT_CHANGED")
    return dict(report=ref,stages=stages,producer=base,validator=validator,
                justification="Same tested producer, pinned environment, historical input and complete scientific gate evidence.")


def prepare(config_path,baseline,baseline_models,destination,validator,*,expected=None,reused=None):
    config=json.loads(Path(config_path).read_text());destination=Path(destination).resolve()
    gate_ref,gate=gate_reference(config)
    baseline_path=verify(baseline);verify(baseline_models)
    historical=json.loads(baseline_path.read_text())
    inventory=set(operator_registry()) if expected is None else set(expected)
    if set(historical["inputs"])!=inventory or not historical["passed"]:
        raise ValueError("EXECUTION_BASELINE_INVENTORY_MISMATCH")
    reused={"1","D3","D5b","D20","71p"} if reused is None else set(reused)
    if not reused<=inventory:raise ValueError("EXECUTION_REUSE_OUTSIDE_INVENTORY")
    producer=config["source_commit"];pilot=Path(config["root"])
    entries={};reuse={};hashes={}
    for operator in sorted(inventory):
        path=Path(config["historical"])/f"op_{operator}.dat"
        hashes[operator]=historical["inputs"][operator]["sha256"]
        if digest_file(path)!=hashes[operator]:raise ValueError("EXECUTION_HISTORICAL_INPUT_CHANGED")
        entries[operator]={};reuse[operator]={}
        for seed in (0,1):
            if operator in reused:
                root=(pilot if operator=="71p" else pilot.parent/"controls"/producer)
                report_path=root/f"seed-{seed}"/operator/"report.json"
                proof=verify_reuse(report_path,operator,seed,producer,hashes[operator],validator,gate)
                reuse[operator][str(seed)]=proof
                entries[operator][str(seed)]=dict(action="reuse",report=proof["report"],stages=STAGES)
            else:
                entries[operator][str(seed)]=dict(action="recompute",stages=STAGES,
                    report_path=str(destination/f"seed-{seed}"/operator/"report.json"),
                    reason="No compatible complete stage data established in the historical reuse inventory.")
    production=dict(config,campaign_id="census-"+producer[:7],root=str(destination),
        queue=str(destination/"queue.sqlite3"),phase="production",seeds=[0,1],
        operators=sorted(inventory-reused),historical_hashes=hashes,worker_limit=4)
    execution=dict(schema_version=1,producer_commit=producer,validator=validator,gate=gate_ref,
        scalar_fermion_only=True,operator_count=len(inventory),seeds=[0,1],entries=entries,
        historical_hashes=hashes,baseline=baseline,baseline_models=baseline_models,
        budget=dict(total_core_hours=3000,includes_validation=True,ledger=config["ledger"]),
        production_config=production)
    # Once published, a repeated preparation must reproduce the identical plan.
    path=destination/"execution-manifest.json"
    with CheckpointStore(destination/"execution-preparation",execution) as store:
        if store.state["complete"]:
            verify(store.state["progress"]["manifest"])
            verify(store.state["progress"]["reuse"])
            verify(store.state["progress"]["config"])
            return path
        atomic_json(destination/"reuse-manifest.json",dict(schema_version=1,gate=gate_ref,operators=reuse))
        atomic_json(destination/"production-config.json",production)
        atomic_json(path,execution)
        store.save(dict(manifest=reference(path),reuse=reference(destination/"reuse-manifest.json"),
                        config=reference(destination/"production-config.json")),complete=True)
    return path


def seal(execution_path,output):
    execution=json.loads(Path(execution_path).read_text())
    config=execution["production_config"];gate_reference(config)
    queue=WorkQueue(config["queue"]);summary=queue.summary()
    if summary["blocked"] or any(row["status"]!="complete" for row in summary["counts"]):
        raise ValueError("CENSUS_EXECUTION_QUEUE_INCOMPLETE")
    entries={}
    for operator,seeds in execution["entries"].items():
        entries[operator]={}
        for seed,entry in seeds.items():
            ref=entry["report"] if entry["action"]=="reuse" else reference(entry["report_path"])
            if entry["action"]=="recompute":
                with queue.connection() as db:
                    unit=db.execute("SELECT status,result FROM units WHERE id=?",(f"{operator}/seed-{seed}/finalise",)).fetchone()
                if (unit is None or unit["status"]!="complete"
                        or json.loads(unit["result"])["report_sha256"]!=ref["sha256"]):
                    raise ValueError("CENSUS_REPORT_NOT_PUBLISHED_BY_OWNING_QUEUE")
            report=json.loads(verify(ref).read_text())
            if (report["source_commit"]!=execution["producer_commit"] or report["operator"]!=operator
                    or str(report["hash_seed"])!=seed
                    or report["historical"]["artifact"]["sha256"]!=execution["historical_hashes"][operator]):
                raise ValueError("CENSUS_REPORT_PRODUCER_OR_INPUT_CHANGED")
            entries[operator][seed]=ref
    result=dict(schema_version=1,producer_commit=execution["producer_commit"],entries=entries,
        execution_manifest=reference(execution_path),baseline=execution["baseline"],baseline_models=execution["baseline_models"])
    # Finalisation's verification stage independently checks every artifact and
    # proves exact cross-seed agreement; sealing itself does not certify totals.
    atomic_json(output,result)
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode",choices=("prepare","seal"));parser.add_argument("input",type=Path)
    parser.add_argument("output",type=Path);parser.add_argument("--baseline",type=Path)
    parser.add_argument("--baseline-models",type=Path)
    args=parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):parser.error("Census verification must run in Slurm")
    if args.mode=="prepare":
        if args.baseline is None or args.baseline_models is None:parser.error("Preparation requires both baseline artifacts")
        provenance=runtime_provenance()
        atomic_json(args.output/"validation-attempts"/(os.environ["SLURM_JOB_ID"]+".json"),provenance)
        result=prepare(args.input,reference(args.baseline),reference(args.baseline_models),args.output,
            dict(environment=scientific_contract(provenance),imported_package=provenance["imported_package"],
                 interpreter=provenance["interpreter"],script_sha256=digest_file(__file__)))
        print(result)
    else:seal(args.input,args.output)
