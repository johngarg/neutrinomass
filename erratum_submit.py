#!/usr/bin/env python3
"""One bounded Slurm scheduling pass; no scientific work runs on the login node."""

import argparse
import fcntl
import json
from pathlib import Path
import shlex
import subprocess
import time

from neutrinomass.database.campaign import CampaignLedger,CampaignBudgetError
from neutrinomass.database.checkpoints import atomic_json
from neutrinomass.database.work_queue import WorkQueue

ACTIVE={"PENDING","RUNNING","CONFIGURING","COMPLETING","SUSPENDED","RESIZING","REQUEUED"}
INFRASTRUCTURE={"NODE_FAIL","PREEMPTED","BOOT_FAIL"}


def accounting(job_ids):
    if not job_ids:
        return {}
    output=subprocess.check_output(["sacct","-n","-P","-j",",".join(job_ids),
        "--format=JobIDRaw,State,AllocCPUS,ElapsedRaw,ExitCode"],text=True)
    result={}
    for line in output.splitlines():
        identifier,status,cpus,elapsed,code=line.split("|")[:5]
        if identifier in job_ids:
            result[identifier]=dict(status=status.split()[0].rstrip("+"),
                cpus=int(cpus),elapsed=int(elapsed),exit_code=code)
    return result


def reconcile(ledger,queue,config_path):
    state=ledger.summary()
    jobs={key:job for key,job in state["jobs"].items()
          if job.get("queue_config")==str(config_path) and "allocated_core_seconds" not in job}
    if any(not job.get("job_id") for job in jobs.values()):
        raise RuntimeError("AMBIGUOUS_SUBMISSION: resolve the reserved Slurm job before submitting again")
    observed=accounting([job["job_id"] for job in jobs.values()])
    for work_id,job in jobs.items():
        record=observed.get(job["job_id"])
        if record is None or record["status"] in ACTIVE:
            continue
        status=record["status"]
        failure=None
        with queue.connection() as db:
            unfinished=[dict(row) for row in db.execute(
                "SELECT id FROM units WHERE worker=? AND status='running'",(job["job_id"],))]
        if status in INFRASTRUCTURE:
            queue.recover_infrastructure(job["job_id"],status)
            failure="infrastructure"
        elif record["exit_code"]=="75:0" and not unfinished:
            status="CONTINUED"
        elif status!="COMPLETED" or record["exit_code"]!="0:0" or unfinished:
            failure="scientific"
            for unit in unfinished:
                queue.fail(unit["id"],job["job_id"],dict(code="UNEXPECTED_SLURM_EXIT",**record))
        ledger.finished(work_id,elapsed_seconds=record["elapsed"],allocated_cpus=record["cpus"],
            status=status,failure_class=failure)
        if failure=="scientific" and not unfinished:
            # A worker may fail before it can claim a unit, e.g. import/pins.
            with ledger.transaction() as state:
                state["submission_block"]=dict(work_id=work_id,record=record)


def launcher(config,config_path,command):
    path=Path(config["root"])/f"launch-{command}.sh"
    path.parent.mkdir(parents=True,exist_ok=True)
    source=shlex.quote(config["source"])
    interpreter=shlex.quote(config["interpreter"])
    config_arg=shlex.quote(str(config_path))
    driver=config.get("worker_driver")
    entry=("-c "+shlex.quote("import runpy; runpy.run_path("+repr(driver)+", run_name='__main__')")
           if driver else "erratum_worker.py")
    path.write_text(f"""#!/bin/bash
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
cd {source}
export PYTHONPATH={source}
export MPLCONFIGDIR=${{SLURM_TMPDIR:-/tmp}}/erratum-mpl-${{SLURM_JOB_ID}}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
{interpreter} -c 'import json,os; p=os.environ.get("SLURM_TMPDIR",os.environ.get("TMPDIR","/tmp")); s=os.statvfs(p); print(json.dumps(dict(event="scratch_capacity",root=p,total_bytes=s.f_blocks*s.f_frsize,free_bytes=s.f_bavail*s.f_frsize)),flush=True)'
exec {interpreter} -u {entry} {command} {config_arg}
""")
    return path


def submit(ledger,config,config_path,*,command,seed):
    with ledger.transaction() as state:
        if state.get("submission_block"):
            raise RuntimeError("SUBMISSION_BLOCKED: diagnose the previous worker failure")
        ordinal=len(state["jobs"])
        retry_of=next((key for key,job in reversed(list(state["jobs"].items()))
            if job.get("queue_config")==str(config_path) and job.get("seed")==seed
            and job.get("failure_class")=="infrastructure" and not job.get("retry_assigned")),None)
    work_id=f"{config['campaign_id']}/{command}/seed-{seed}/{ordinal:06d}"
    ledger.reserve(work_id,phase=config["phase"],retry_of=retry_of)
    with ledger.transaction() as state:
        state["jobs"][work_id].update(queue_config=str(config_path),seed=seed,command=command)
        if retry_of is not None:
            state["jobs"][retry_of]["retry_assigned"]=work_id
    script=launcher(config,config_path,command)
    logs=Path(config["root"])/"logs"
    logs.mkdir(parents=True,exist_ok=True)
    process=subprocess.run(["sbatch","--parsable",f"--comment={work_id}",
        f"--job-name=erratum-{seed}",f"--chdir={config['source']}",
        f"--output={logs}/%j.log",f"--export=ALL,PYTHONHASHSEED={seed}",str(script)],
        text=True,capture_output=True)
    if process.returncode:
        ledger.finished(work_id,elapsed_seconds=0,allocated_cpus=1,status="SUBMISSION_FAILED",failure_class="scientific")
        with ledger.transaction() as state:
            state["submission_block"]=dict(work_id=work_id,stderr=process.stderr)
        raise RuntimeError(f"Slurm submission failed: {process.stderr}")
    job_id=process.stdout.strip().split(";")[0]
    if not job_id.isdigit():
        raise RuntimeError("AMBIGUOUS_SUBMISSION: Slurm did not return a job identifier")
    ledger.submitted(work_id,job_id)
    return dict(work_id=work_id,job_id=job_id,seed=seed)


def _tick(config_path,*,prepare=False):
    config_path=Path(config_path).resolve()
    config=json.loads(config_path.read_text())
    ledger=CampaignLedger(config["ledger"])
    queue=WorkQueue(config["queue"])
    reconcile(ledger,queue,config_path)
    summary=queue.summary()
    if summary["blocked"]:
        return dict(status="blocked",queue=summary,accounting=ledger.summary())
    if prepare:
        already=any(job.get("queue_config")==str(config_path) and job.get("command")=="prepare"
                    for job in ledger.summary()["jobs"].values())
        if already or summary["counts"] or (Path(config["root"])/"prepared.json").exists():
            raise ValueError("Preparation already started; inspect or use tick")
        return submit(ledger,config,config_path,command="prepare",seed=config["seeds"][0])
    if not (Path(config["root"])/"prepared.json").exists():
        return dict(status="awaiting-preparation",accounting=ledger.summary())
    remaining=[row for row in summary["counts"] if row["status"]!="complete"]
    if not remaining:
        return dict(status="complete",queue=summary,accounting=ledger.summary())
    # Complete seed zero independently before starting seed one.
    seed=min(row["seed"] for row in remaining)
    submitted=[]
    while True:
        with queue.connection() as db:
            ready=db.execute("""SELECT COUNT(*) FROM units u WHERE status='pending' AND seed=?
                AND NOT EXISTS (SELECT 1 FROM dependencies d LEFT JOIN units p ON p.id=d.dependency
                WHERE d.unit=u.id AND (p.status IS NULL OR p.status!='complete'))""",(seed,)).fetchone()[0]
        live=[j for j in ledger.summary()["jobs"].values() if "allocated_core_seconds" not in j]
        same=sum(j.get("queue_config")==str(config_path) and j.get("seed")==seed for j in live)
        if not ready or same>=ready or len(live)>=4:
            break
        try:
            submitted.append(submit(ledger,config,config_path,command="work",seed=seed))
        except CampaignBudgetError as error:
            return dict(status="budget-stop",reason=str(error),submitted=submitted,accounting=ledger.summary())
    return dict(status="running",submitted=submitted,queue=summary,accounting=ledger.summary())


def tick(config_path,*,prepare=False):
    config_path=Path(config_path).resolve()
    with config_path.with_suffix(".controller.lock").open("a+") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        return _tick(config_path,prepare=prepare)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config",type=Path)
    parser.add_argument("--prepare",action="store_true")
    args=parser.parse_args()
    print(json.dumps(tick(args.config,prepare=args.prepare),sort_keys=True,indent=2))
