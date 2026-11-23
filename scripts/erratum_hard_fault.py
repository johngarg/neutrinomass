"""Kill a replica after a durable checkpoint, then recover it in a new process."""
import argparse
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

from neutrinomass.database.checkpoints import atomic_json,digest_file


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind",choices=("generation","coverage"))
    parser.add_argument("config",type=Path);parser.add_argument("reference",type=Path)
    parser.add_argument("replica",type=Path);parser.add_argument("--unit",type=int,required=True)
    args=parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Fault acceptance must run in Slurm")
    script=Path(__file__).with_name("erratum_recovery.py")
    command=[sys.executable,"-u",str(script),args.kind,str(args.config),str(args.reference),
             str(args.replica),"--unit",str(args.unit)]
    checkpoint=(args.replica/"units"/f"{args.unit:08d}"/"generation/checkpoint.json"
                if args.kind=="generation" else args.replica/"checkpoint.json")
    archive=args.replica/"interruptions"/os.environ["SLURM_JOB_ID"]
    evidence=archive.with_suffix(".json")
    previous=[]
    for path in (args.replica/"interruptions").glob("*.json"):
        item=json.loads(path.read_text())
        if (item.get("method")=="SIGKILL" and item.get("injected")
                and item.get("kind")==args.kind and item.get("unit")==args.unit):
            for saved in item["checkpoints"]:
                if digest_file(saved["path"])!=saved["sha256"]:
                    raise ValueError("HARD_FAULT_EVIDENCE_CORRUPT")
            previous.append(path)
    child=None;stopping=False
    def stop(signum,frame):
        nonlocal stopping
        stopping=True
        if child is not None and child.poll() is None:child.send_signal(signal.SIGUSR1)
    signal.signal(signal.SIGUSR1,stop);signal.signal(signal.SIGTERM,stop)
    if not previous:
        if checkpoint.exists():
            parser.error("Use a fresh replica for a new hard-fault experiment")
        child=subprocess.Popen(command+["--checkpoint-seconds","1"])
        started=time.monotonic()
        while child.poll() is None:
            if stopping:
                child.wait();return 75
            if checkpoint.exists():
                data=json.loads(checkpoint.read_text())
                if not data["complete"] and data["progress"]:
                    # Let the child perform work after this checkpoint; that
                    # uncommitted work is deliberately lost by SIGKILL.
                    time.sleep(0.2)
                    if child.poll() is None:
                        child.kill();child.wait()
                        data=json.loads(checkpoint.read_text())
                        archive.mkdir(parents=True,exist_ok=True)
                        shutil.copyfile(checkpoint,archive/"checkpoint.json")
                        for artifact in data["artifacts"].values():
                            shutil.copyfile(checkpoint.parent/artifact["path"],archive/artifact["path"])
                        atomic_json(evidence,dict(kind=args.kind,unit=args.unit,injected=True,
                            method="SIGKILL",child_returncode=child.returncode,
                            job_id=os.environ["SLURM_JOB_ID"],validator_sha256=digest_file(__file__),
                            checkpoints=[dict(path=str(archive/"checkpoint.json"),sha256=digest_file(archive/"checkpoint.json"),
                                              progress=data["progress"])]))
                        break
            if time.monotonic()-started>300:
                child.send_signal(signal.SIGUSR1);child.wait()
                raise RuntimeError("Replica did not publish a checkpoint within five minutes")
            time.sleep(0.1)
        if not evidence.exists():
            raise RuntimeError(f"Replica exited before a hard fault could be injected: {child.returncode}")
    # A fresh interpreter resumes only the immutable chunks and sealed SQLite
    # snapshot. No Python or node-local transaction state survives the kill.
    if stopping:return 75
    child=subprocess.Popen(command)
    return child.wait()


if __name__=="__main__":
    raise SystemExit(main())
