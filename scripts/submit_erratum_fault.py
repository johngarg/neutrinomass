"""Reserve one bounded Slurm slot for a pilot fault replica in the active task."""
import argparse
import fcntl
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import time

from neutrinomass.database.campaign import CampaignLedger,CampaignBudgetError
from neutrinomass.database.checkpoints import atomic_json,digest_file


def worker_limit(path,limit):
    config=json.loads(path.read_text())
    digest=digest_file(path)
    archive=path.parent/"config-revisions"/(digest+".json")
    archive.parent.mkdir(exist_ok=True)
    if not archive.exists():shutil.copyfile(path,archive)
    if digest_file(archive)!=digest:
        raise ValueError("CONFIG_ARCHIVE_HASH_MISMATCH")
    config["worker_limit"]=limit
    atomic_json(path,config)
    return config


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config",type=Path)
    parser.add_argument("kind",choices=("generation","reduction","coverage"))
    parser.add_argument("--unit",type=int,required=True)
    parser.add_argument("--seed",type=int,choices=(0,1),default=0)
    args=parser.parse_args()
    path=args.config.resolve()
    config=worker_limit(path,3)
    ledger=CampaignLedger(config["ledger"])
    root=Path(config["root"])
    replica=root/"recovery"/f"{args.kind}-seed-{args.seed}-unit-{args.unit}"
    work_id=f"{config['campaign_id']}/hard-fault/{args.kind}/seed-{args.seed}/unit-{args.unit}"
    reference=root/f"seed-{args.seed}"/"71p"
    checkpoint=(reference/"units"/f"{args.unit:08d}"/"structural/checkpoint.json"
                if args.kind=="generation" else reference/"structural/checkpoint.json" if args.kind=="reduction"
                else reference/"coverage/shards"/f"{args.unit:08d}"/"checkpoint.json")
    if not json.loads(checkpoint.read_text())["complete"]:
        raise ValueError("REFERENCE_UNIT_NOT_COMPLETE")
    print(json.dumps(dict(event="waiting_for_fault_slot",work_id=work_id)),flush=True)
    controller_lock=path.with_suffix(".controller.lock").open("a+")
    while True:
        try:
            fcntl.flock(controller_lock,fcntl.LOCK_EX)
            ledger.reserve(work_id,phase="validation")
            break
        except CampaignBudgetError as error:
            fcntl.flock(controller_lock,fcntl.LOCK_UN)
            if str(error)!="Four-worker concurrency limit":raise
            time.sleep(10)
    with ledger.transaction() as state:
        state["jobs"][work_id].update(queue_config=str(path),seed=args.seed,
            command="recovery-"+args.kind,replica=str(replica),validator_sha256=digest_file(__file__))
    q=shlex.quote
    driver=Path(__file__).resolve().with_name("erratum_hard_fault.py")
    launcher=root/f"launch-hard-fault-{args.kind}-{args.seed}-{args.unit}.sh"
    launcher.write_text(f"""#!/bin/bash
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
cd {q(config['source'])}
export PYTHONPATH={q(config['source'])}
export PYTHONHASHSEED={args.seed}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR=${{SLURM_TMPDIR:-/tmp}}/erratum-fault-mpl-${{SLURM_JOB_ID}}
exec {q(config['interpreter'])} -u {q(str(driver))} {args.kind} {q(str(path))} {q(str(reference))} {q(str(replica))} --unit {args.unit}
""")
    result=subprocess.run(["sbatch","--parsable",f"--comment={work_id}",
        "--job-name=erratum-hard-fault",f"--output={root}/logs/%j.log",str(launcher)],
        text=True,capture_output=True)
    if result.returncode:
        ledger.finished(work_id,elapsed_seconds=0,allocated_cpus=1,status="SUBMISSION_FAILED",failure_class="scientific")
        with ledger.transaction() as state:
            state["submission_block"]=dict(work_id=work_id,stderr=result.stderr)
        raise RuntimeError(result.stderr)
    job_id=result.stdout.strip().split(";")[0]
    if not job_id.isdigit():raise RuntimeError("AMBIGUOUS_FAULT_SUBMISSION")
    ledger.submitted(work_id,job_id)
    worker_limit(path,4)
    fcntl.flock(controller_lock,fcntl.LOCK_UN)
    controller_lock.close()
    print(json.dumps(dict(event="fault_submitted",work_id=work_id,job_id=job_id,replica=str(replica))),flush=True)


if __name__=="__main__":main()
