"""Compare global model unions without mistaking operator transfers for changes.

Run after complete census publication. The published identification convention
is used for both unfiltered and historically filtered representation sets.
Species-aware mechanisms and exact Lagrangians retain their separate reports.
"""
import argparse
import ast
import json
import os
from pathlib import Path

from neutrinomass.database.checkpoints import CheckpointStore,WorkBudget,BudgetExpired,atomic_json,canonical_json,digest_file
from neutrinomass.database.filtering import published_model_key
from neutrinomass.database.provenance import runtime_provenance,scientific_contract
from neutrinomass.database.resumable_reduction import lines_from_paths


def verified(ref):
    path=Path(ref["path"])
    if not path.is_file() or digest_file(path)!=ref["sha256"]:raise ValueError("GLOBAL_COMPARISON_INPUT_CHANGED")
    return path


def sources_from_report(report_path,input_manifest):
    report=json.loads(Path(report_path).read_text());manifest=json.loads(Path(input_manifest).read_text())
    if not report["complete_census"] or report["operator_count"]!=243 or report["seeds"]!=[0,1]:
        raise ValueError("GLOBAL_COMPARISON_REQUIRES_COMPLETE_CENSUS")
    if report["source_reports"]!=manifest["entries"]:raise ValueError("GLOBAL_COMPARISON_MANIFEST_CHANGED")
    sources=[]
    for operator,seeds in sorted(report["source_reports"].items()):
        current=json.loads(verified(seeds["0"]).read_text())
        sources.extend([
            dict(mode="unfiltered",side="historical",kind="legacy",operator=operator,artifact=current["historical"]["artifact"]),
            dict(mode="unfiltered",side="current",kind="json",operator=operator,artifact=current["artifacts"]["democratic_models"])])
    checkpoint=verified(report["dependencies"]["published"])
    stage=json.loads(checkpoint.read_text())
    with CheckpointStore(checkpoint.parent,stage["contract"]) as store:
        if not store.state["complete"]:raise ValueError("GLOBAL_COMPARISON_FILTER_INCOMPLETE")
        artifact=store.state["artifacts"]["output"]
        published=dict(path=str(checkpoint.parent/artifact["path"]),sha256=artifact["sha256"])
    sources.extend([
        dict(mode="published",side="historical",kind="json",operator=None,artifact=manifest["baseline_models"]),
        dict(mode="published",side="current",kind="json",operator=None,artifact=published)])
    return sources


def compare(sources,directory,contract,*,budget=None,scratch=None):
    budget=budget or WorkBudget();directory=Path(directory)
    for source in sources:verified(source["artifact"])
    identity=dict(contract,algorithm="global-model-union-comparison-v1",sources=sources,
                  identification="published",scale_unit="TeV")
    with CheckpointStore(directory,identity) as store:
        if store.state["complete"]:return store.state
        with store.database(scratch) as db:
            db.execute("CREATE TABLE IF NOT EXISTS models(mode TEXT,side TEXT,key TEXT,fields TEXT,PRIMARY KEY(mode,side,key))")
            progress=dict(store.state["progress"] or dict(source_index=0,record=0,scanned=0))
            def save():store.snapshot(db,progress)
            for index in range(progress["source_index"],len(sources)):
                source=sources[index]
                for ordinal,line in enumerate(lines_from_paths([source["artifact"]["path"]]),1):
                    if ordinal<=progress["record"]:continue
                    if source["kind"]=="legacy":
                        prefix="LazyCompletion(head="
                        if not line.startswith(prefix):raise ValueError("GLOBAL_COMPARISON_LEGACY_HEADER")
                        head=ast.literal_eval(line[len(prefix):].split(", tail=",1)[0])
                        if head["operator_name"]!=source["operator"]:raise ValueError("GLOBAL_COMPARISON_OPERATOR_CHANGED")
                        fields=tuple(sorted(set(head["quantum_numbers"])))
                    else:fields=tuple(json.loads(line)["fields"])
                    key=canonical_json(published_model_key(fields))
                    db.execute("INSERT OR IGNORE INTO models VALUES (?,?,?,?)",
                        (source["mode"],source["side"],key,canonical_json(fields)))
                    progress.update(record=ordinal,scanned=progress["scanned"]+1)
                    budget.boundary(save,progress=dict(stage="global-model-comparison",**progress))
                progress.update(source_index=index+1,record=0)
            counts={}
            for mode in ("unfiltered","published"):
                sizes={side:db.execute("SELECT COUNT(*) FROM models WHERE mode=? AND side=?",(mode,side)).fetchone()[0]
                       for side in ("historical","current")}
                common=db.execute("SELECT COUNT(*) FROM models a JOIN models b ON a.mode=b.mode AND a.key=b.key "
                                  "WHERE a.mode=? AND a.side='historical' AND b.side='current'",(mode,)).fetchone()[0]
                counts[mode]=dict(sizes,unchanged=common,added=sizes["current"]-common,removed=sizes["historical"]-common)
            output=store.publish_lines(canonical_json(dict(mode=mode,code="MODEL_ADDED" if side=="current" else "MODEL_REMOVED",
                key=json.loads(key),fields=json.loads(fields))) for mode,side,key,fields in db.execute(
                    "SELECT a.mode,a.side,a.key,a.fields FROM models a LEFT JOIN models b "
                    "ON a.mode=b.mode AND a.key=b.key AND a.side!=b.side WHERE b.key IS NULL ORDER BY a.mode,a.side,a.key"))
            progress["modes"]=counts
            store.save(progress,artifacts={"output":output})
            store.snapshot(db,progress,complete=True)
            return store.state


def validate_totals(counts,report,baseline):
    expected={
        "unfiltered":dict(historical=baseline["conjugacy_unfiltered_models"],
                          current=report["global_representation_sets"]),
        "published":dict(historical=baseline["filtered_models"],
                         current=report["filtered_published"]["representation_sets"])}
    for mode,sides in expected.items():
        for side,total in sides.items():
            if counts[mode][side]!=total:
                raise ValueError("GLOBAL_COMPARISON_TOTAL_MISMATCH: "+mode+"/"+side)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report",type=Path);parser.add_argument("input_manifest",type=Path);parser.add_argument("output",type=Path)
    args=parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):parser.error("Global census comparison must run in Slurm")
    p=runtime_provenance();contract=dict(scientific_contract(p),validator_sha256=digest_file(__file__),
        report_sha256=digest_file(args.report),manifest_sha256=digest_file(args.input_manifest))
    budget=WorkBudget(3300);budget.scratch_root=os.environ.get("SLURM_TMPDIR")
    try:
        with budget.signals():
            state=compare(sources_from_report(args.report,args.input_manifest),args.output,contract,
                          budget=budget,scratch=budget.scratch_root)
    except BudgetExpired:raise SystemExit(75)
    manifest=json.loads(args.input_manifest.read_text())
    baseline=json.loads(verified(manifest["baseline"]).read_text())
    validate_totals(state["progress"]["modes"],json.loads(args.report.read_text()),baseline)
    result=dict(schema_version=2,scale_unit="TeV",complete=True,identification="published",
        global_model_changes=state["progress"]["modes"],producer=p,
        source_report=dict(path=str(args.report),sha256=contract["report_sha256"]),
        input_manifest=dict(path=str(args.input_manifest),sha256=contract["manifest_sha256"]),
        historical_baseline=manifest["baseline"],
        differences=dict(path=str(args.output/state["artifacts"]["output"]["path"]),sha256=state["artifacts"]["output"]["sha256"]),
        checkpoint=dict(path=str(args.output/"checkpoint.json"),sha256=digest_file(args.output/"checkpoint.json")),
        low_scale_genuineness_proved=False)
    atomic_json(args.output/"report.json",result)
    print(json.dumps(result,sort_keys=True))
