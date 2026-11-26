"""Follow the active pilot, recovery gate, census and finalisation to completion.

Only metadata, accounting and Slurm submission run on the login node. Every
scientific command runs in a bounded allocation. Restarting this driver reads
the same ledger and queues; it cannot duplicate an ambiguous submission.
"""
import argparse
import fcntl
import json
import os
from pathlib import Path
import shlex
import subprocess
import time

from erratum_submit import tick,reconcile
from neutrinomass.database.campaign import CampaignLedger,CampaignBudgetError
from neutrinomass.database.checkpoints import atomic_json,digest_file,identity_hash
from neutrinomass.database.work_queue import WorkQueue,progress_signature
from scripts.submit_erratum_fault import worker_limit


def verified(ref):
    path=Path(ref["path"])
    if not path.is_file() or digest_file(path)!=ref["sha256"]:
        raise ValueError("CHAIN_REFERENCE_CHANGED: "+str(path))
    return path


def complete(path):
    return Path(path).is_file() and json.loads(Path(path).read_text()).get("complete") is True


def external(ledger,pilot_path,chain,stage,spec):
    """Reserve and submit one immutable command, or resume its saved work."""
    pilot=json.loads(pilot_path.read_text())
    identity=identity_hash(spec)
    jobs=[(key,job) for key,job in ledger.summary()["jobs"].items()
          if job.get("chain_stage")==stage and job.get("chain_id")==chain["id"]]
    if any(job["chain_identity"]!=identity for _,job in jobs):
        raise ValueError("CHAIN_COMMAND_CHANGED_REQUIRES_DIAGNOSIS")
    if jobs:
        old_id,old=jobs[-1]
        if "allocated_core_seconds" not in old:return "running"
        if old["status"]=="COMPLETED":
            for path in spec["evidence"]:
                if not Path(path).is_file():raise ValueError("CHAIN_SUCCESS_WITHOUT_EVIDENCE")
            with ledger.transaction() as state:
                job=state["jobs"][old_id]
                if "evidence_references" not in job:
                    job["evidence_references"]=[dict(path=str(path),sha256=digest_file(path)) for path in spec["evidence"]]
                for ref in job["evidence_references"]:verified(ref)
            return "complete"
        if old.get("failure_class")=="infrastructure":
            retry_of=old_id
        elif old["status"]=="CONTINUED":
            retry_of=None
            signature=identity_hash({str(path):progress_signature(path) for path in spec["progress"]})
            with ledger.transaction() as state:
                job=state["jobs"][old_id]
                if "continuation_signature" not in job:
                    previous=jobs[-2][1] if len(jobs)>1 else {}
                    repeated=previous.get("no_progress",0)+1 if previous.get("continuation_signature")==signature else 0
                    job.update(continuation_signature=signature,no_progress=repeated)
                if job["no_progress"]>=2:
                    state["submission_block"]=dict(work_id=old_id,code="CHAIN_REPEATED_NO_PROGRESS")
            if ledger.summary().get("submission_block"):raise ValueError("CHAIN_REPEATED_NO_PROGRESS")
        else:raise ValueError("CHAIN_FAILED_UNIT_REQUIRES_DIAGNOSIS")
    else:retry_of=None
    work_id=f"{chain['id']}/{stage}/{len(jobs):04d}"
    # Use the same lock as the ordinary pilot controller. The reservation and
    # Slurm job binding form one short serialized submission operation.
    with pilot_path.with_suffix(".controller.lock").open("a+") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if ledger.summary().get("submission_block"):raise ValueError("CHAIN_SUBMISSION_BLOCKED")
        try:ledger.reserve(work_id,phase=spec["phase"],retry_of=retry_of)
        except CampaignBudgetError as error:
            if str(error)=="Four-worker concurrency limit":return "awaiting-slot"
            raise
        with ledger.transaction() as state:
            state["jobs"][work_id].update(queue_config=str(pilot_path),seed=0,command="chain-"+stage,
                chain_id=chain["id"],chain_stage=stage,chain_identity=identity,specification=spec)
        root=Path(chain["root"]);root.mkdir(parents=True,exist_ok=True)
        batch=root/(stage+f"-{len(jobs):04d}.sh");q=shlex.quote
        command=" ".join(q(str(x)) for x in (pilot["interpreter"],"-u",*spec["argv"]))
        batch.write_text(f"""#!/bin/bash
#SBATCH --partition=sapphire
#SBATCH --account=punim0011
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --hint=nomultithread
#SBATCH --signal=B:USR1@1800
set -euo pipefail
module purge
module load GCCcore/11.3.0 Python/3.11.3
cd {q(spec['source'])}
export PYTHONPATH={q(spec['source'])} PYTHONHASHSEED=0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR=${{SLURM_TMPDIR:-/tmp}}/erratum-chain-${{SLURM_JOB_ID}}
exec {command}
""")
        result=subprocess.run(["sbatch","--parsable",f"--comment={work_id}",
            "--job-name=erratum-"+stage,f"--output={root}/%j.log",str(batch)],text=True,capture_output=True)
        if result.returncode:
            ledger.finished(work_id,elapsed_seconds=0,allocated_cpus=1,status="SUBMISSION_FAILED",failure_class="scientific")
            with ledger.transaction() as state:state["submission_block"]=dict(work_id=work_id,stderr=result.stderr)
            raise RuntimeError(result.stderr)
        job_id=result.stdout.strip().split(";")[0]
        if not job_id.isdigit():raise ValueError("CHAIN_AMBIGUOUS_SUBMISSION")
        ledger.submitted(work_id,job_id)
        print(json.dumps(dict(event="chain_submitted",stage=stage,job_id=job_id,work_id=work_id)),flush=True)
    return "running"


def queue_tick(path):
    config=json.loads(path.read_text())
    queue=WorkQueue(config["queue"])
    state=CampaignLedger(config["ledger"]).summary()
    prepared=any(job.get("queue_config")==str(path) and job.get("command")=="prepare" for job in state["jobs"].values())
    if not prepared and not queue.summary()["counts"] and not (Path(config["root"])/"prepared.json").exists():
        return dict(status="preparing",submission=tick(path,prepare=True))
    return tick(path)


def step(chain):
    for ref in chain["orchestrator_sources"]:verified(ref)
    pilot_path=Path(chain["pilot_config"]).resolve();pilot=json.loads(pilot_path.read_text())
    ledger=CampaignLedger(pilot["ledger"]);queue=WorkQueue(pilot["queue"])
    reconcile(ledger,queue,pilot_path)
    if ledger.summary().get("submission_block"):raise ValueError("CHAIN_SUBMISSION_BLOCKED")
    if queue.summary()["blocked"]:raise ValueError("CHAIN_PILOT_SCIENTIFIC_FAILURE")
    scripts=Path(chain["scripts"]);root=Path(pilot["root"])
    reference=root/"seed-0/71p"
    # These are exact certificates from the independently completed full tests.
    for ref in chain["finaliser_validation"]:
        data=json.loads(verified(ref).read_text())
        if (not data["full_suite_passed"] or data["environment"]["source_commit"]!=chain["finaliser_commit"]
                or data["environment"]["dirty_worktree"]):
            raise ValueError("CHAIN_FINALISER_VALIDATION_FAILED")
    if {json.loads(verified(ref).read_text())["environment"]["hash_seed"] for ref in chain["finaliser_validation"]}!={"0","1"}:
        raise ValueError("CHAIN_FINALISER_VALIDATION_SEEDS_INCOMPLETE")
    generation=Path(chain["generation_recovery"])
    recovered=json.loads((generation/"recovery-evidence.json").read_text())
    if not recovered["complete"] or not recovered["ordered_exact_agreement"]:
        raise ValueError("CHAIN_GENERATION_RECOVERY_FAILED")
    recoveries=[generation]
    pending=[]
    for kind in ("reduction","coverage"):
        replica=root/"recovery"/f"{kind}-seed-0-unit-0";recoveries.append(replica)
        ready=(complete(reference/"structural/checkpoint.json") if kind=="reduction" else (reference/"report.json").is_file())
        pending.append((kind,ready,dict(phase="validation",source=pilot["source"],
            argv=[str(scripts/"erratum_hard_fault.py"),kind,str(pilot_path),str(reference),str(replica),"--unit","0"],
            progress=[str(replica)],evidence=[str(replica/"recovery-evidence.json")])))
    replica=root/"recovery/substitution-seed-0-unit-0";recoveries.append(replica)
    pending.append(("substitution",(reference/"report.json").is_file(),dict(phase="validation",source=pilot["source"],
        argv=[str(scripts/"erratum_recovery.py"),"substitution",str(pilot_path),str(reference),str(replica),
              "--unit","0","--recovered-generation",str(generation)],
        progress=[str(replica)],evidence=[str(replica/"recovery-evidence.json")])))
    waiting=False
    for stage,ready,spec in pending:
        if not ready:waiting=True;break
        status=external(ledger,pilot_path,chain,stage,spec)
        if status!="complete":
            worker_limit(pilot_path,3)
            result=tick(pilot_path)
            if result["status"] in {"blocked","budget-stop"}:return result
            return dict(status="pilot-and-recovery",stage=stage,external=status,pilot=result["status"])
    worker_limit(pilot_path,4)
    result=tick(pilot_path)
    if result["status"]!="complete":return result
    if waiting:raise ValueError("CHAIN_COMPLETED_PILOT_LACKS_REQUIRED_STAGE")
    # Reconcile all earlier allocations before certifying pilot resource costs.
    if any("allocated_core_seconds" not in job and job.get("queue_config")==str(pilot_path)
           and job.get("command") in {"prepare","work"} for job in ledger.summary()["jobs"].values()):
        return dict(status="awaiting-accounting")
    acceptance=root/"acceptance"
    gate_script=verified(chain["gate_script"]) if "gate_script" in chain else scripts/"validate_71p_gate.py"
    command=[str(gate_script),str(pilot_path),"--baseline",str(verified(chain["baseline"]))]
    for replica in recoveries:command.extend(["--recovery",str(replica)])
    status=external(ledger,pilot_path,chain,"gate",dict(phase="validation",source=pilot["source"],
        argv=command,progress=[str(acceptance)],evidence=[str(acceptance/"gate.json")]))
    if status!="complete":return dict(status="gate",external=status)
    production=Path(chain["production_root"]);execution=production/"execution-manifest.json"
    status=external(ledger,pilot_path,chain,"execution-plan",dict(phase="production",source=pilot["source"],
        argv=[str(scripts/"prepare_erratum_campaign.py"),"prepare",str(pilot_path),str(production),
              "--baseline",str(verified(chain["baseline"])),"--baseline-models",str(verified(chain["baseline_models"]))],
        progress=[str(production/"execution-preparation")],evidence=[str(execution),str(production/"production-config.json")]))
    if status!="complete":return dict(status="execution-plan",external=status)
    result=queue_tick(production/"production-config.json")
    if result["status"]!="complete":return dict(result,phase="production")
    manifest=production/"finalisation-inputs.json"
    status=external(ledger,pilot_path,chain,"seal-census",dict(phase="finalisation",source=pilot["source"],
        argv=[str(scripts/"prepare_erratum_campaign.py"),"seal",str(execution),str(manifest)],
        progress=[],evidence=[str(manifest)]))
    if status!="complete":return dict(status="seal-census",external=status)
    final=Path(chain["finalisation_root"]);config_path=final/"config.json"
    if not config_path.exists():
        atomic_json(config_path,dict(root=str(final),queue=str(final/"queue.sqlite3"),ledger=pilot["ledger"],
            campaign_id=chain["id"]+"-finalisation",source=chain["finaliser_source"],source_commit=chain["finaliser_commit"],
            worker_driver=str(Path(chain["finaliser_source"])/"erratum_finalise_worker.py"),interpreter=pilot["interpreter"],
            phase="finalisation",seeds=[0],worker_limit=4,manifest=str(manifest),manifest_sha256=digest_file(manifest)))
    result=queue_tick(config_path)
    if result["status"]!="complete":return dict(result,phase="finalisation")
    report=json.loads((final/"report.json").read_text())
    if not report["complete_census"] or report["operator_count"]!=243 or report["seeds"]!=[0,1]:
        raise ValueError("CHAIN_FINAL_CENSUS_INCOMPLETE")
    if any("allocated_core_seconds" not in job for job in ledger.summary()["jobs"].values()):
        return dict(status="awaiting-final-accounting")
    return dict(status="complete",report=dict(path=str(final/"report.json"),sha256=digest_file(final/"report.json")))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("config",type=Path)
    args=parser.parse_args();chain=json.loads(args.config.read_text())
    root=Path(chain["root"]);root.mkdir(parents=True,exist_ok=True)
    with (root/"driver.lock").open("a+") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        previous=None
        while True:
            try:result=step(chain)
            except BlockingIOError:
                time.sleep(1);continue
            except CampaignBudgetError as error:result=dict(status="budget-stop",reason=str(error))
            except Exception as error:
                result=dict(status="blocked",error_type=type(error).__name__,reason=str(error))
            pilot=json.loads(Path(chain["pilot_config"]).read_text())
            accounting=CampaignLedger(pilot["ledger"]).summary()
            snapshot=dict({k:v for k,v in result.items() if k!="accounting"},updated_unix=time.time(),driver_sha256=digest_file(__file__),
                config_sha256=digest_file(args.config),charged_core_hours=accounting["charged_core_hours"],
                validation_charged_core_hours=accounting["validation_charged_core_hours"])
            atomic_json(root/"status.json",snapshot)
            event=(result["status"],result.get("phase"),result.get("stage"))
            if event!=previous:print(json.dumps(snapshot,sort_keys=True),flush=True);previous=event
            if result["status"] in {"complete","blocked","budget-stop"}:return 0 if result["status"]=="complete" else 2
            time.sleep(30)


if __name__=="__main__":raise SystemExit(main())
