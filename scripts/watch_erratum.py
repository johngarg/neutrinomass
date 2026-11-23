"""Monitor this active campaign; submit only budgeted work from its fixed queue."""
import argparse
import json
from pathlib import Path
import time

from erratum_submit import tick
from neutrinomass.database.checkpoints import atomic_json,digest_file


def watch(config_path):
    config_path=Path(config_path).resolve()
    root=Path(json.loads(config_path.read_text())["root"])
    while True:
        result=tick(config_path)
        accounting=result["accounting"]
        snapshot=dict(status=result["status"],submitted=result.get("submitted",[]),
            queue=result.get("queue"),charged_core_hours=accounting["charged_core_hours"],
            validation_charged_core_hours=accounting["validation_charged_core_hours"],
            remaining_core_hours=accounting["remaining_core_hours"],updated_unix=time.time(),
            watcher_sha256=digest_file(__file__))
        if "reason" in result:snapshot["reason"]=result["reason"]
        atomic_json(root/"watch.json",snapshot)
        if snapshot["submitted"] or snapshot["status"]!="running":
            print(json.dumps(snapshot,sort_keys=True),flush=True)
        if snapshot["status"] in {"complete","blocked","budget-stop"}:
            return 0 if snapshot["status"]=="complete" else 2
        time.sleep(30)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config",type=Path)
    args=parser.parse_args()
    raise SystemExit(watch(args.config))
