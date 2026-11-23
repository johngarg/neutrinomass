"""Inventory old stage artifacts without treating reports as completion data."""
import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess

from neutrinomass.database.checkpoints import digest_file,atomic_json


STAGES={"generation":"generator","structural-reduction":"structural_exact_classes",
        "amplitudes":"exact_classes","historical-coverage":"historical_coverage",
        "global-lagrangians":"global_lagrangians"}


def inventory(archives,historical,held,output):
    validator=subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip()
    rows=[]
    for archive in archives:
        for path in sorted(Path(archive).rglob("report.json")):
            report=json.loads(path.read_text())
            for stage,key in STAGES.items():
                entry=report.get("artifacts",{}).get(key)
                row=dict(report=str(path),report_sha256=digest_file(path),
                    operator=report.get("operator"),seed=report.get("hash_seed"),stage=stage,
                    producer=report.get("source_commit"),validator=validator,reuse=False)
                if not isinstance(entry,dict) or not entry.get("path"):
                    row["code"]="REQUIRED_ARTIFACT_NOT_RECORDED"
                else:
                    recorded=Path(entry["path"])
                    candidates=(recorded,path.parent/recorded.name,path.parent/(recorded.name+".gz"))
                    actual=next((p for p in candidates if p.is_file()),None)
                    row["producer_reference"]=entry
                    if actual is None:
                        row["code"]="REQUIRED_ARTIFACT_MISSING"
                    else:
                        actual_hash=digest_file(actual)
                        row["artifact"]=dict(path=str(actual),sha256=actual_hash,bytes=actual.stat().st_size)
                        # Packaging may have compressed the producer's file;
                        # only an explicit matching stored digest certifies it.
                        row["code"]=("SCIENTIFIC_COMPATIBILITY_NOT_ESTABLISHED"
                            if actual_hash==entry.get("sha256") else "PACKAGING_OR_HASH_NOT_VERIFIED")
                rows.append(row)
    historical_rows=[dict(operator=p.stem[3:],path=str(p),sha256=digest_file(p),
        stage="historical-input",reuse=True,code="IMMUTABLE_INPUT_HASH_VERIFIED",validator=validator)
        for p in sorted(Path(historical).glob("op_*.dat"))]
    result=dict(schema_version=1,validator_commit=validator,stages=rows,historical_inputs=historical_rows,
        held_reference=dict(path=str(held),reuse=False,code="FRESH_SPARTAN_71P_REQUIRED"),
        summary=dict(Counter(row["code"] for row in rows)),
        policy="Existing bytes are preserved. Matching hashes establish integrity, not scientific compatibility. No unverified stage is reused.")
    atomic_json(output,result)
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive",type=Path,action="append",required=True)
    parser.add_argument("--historical",type=Path,required=True)
    parser.add_argument("--held",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    report=inventory(args.archive,args.historical,args.held,args.output)
    print(json.dumps(report["summary"],indent=2))
