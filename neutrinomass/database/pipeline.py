"""One staged implementation for individual work units and whole operators."""

from collections import Counter,defaultdict
import hashlib
import json
from pathlib import Path
import shutil
import sqlite3
import tempfile
import time

from neutrinomass.database.checkpoints import CheckpointStore, WorkBudget, digest_file, atomic_json
from neutrinomass.database.planning import generation_plan
from neutrinomass.database.stages import (
    generate_unit,map_records,output_paths,decode_record,amplitude_record,historical_transform,
)
from neutrinomass.database.resumable_reduction import reduce_resumable,lines_from_paths,input_contract
from neutrinomass.database.coverage import historical_coverage
from neutrinomass.database.rebuild import (
    operator_registry,operator_kind,operator_scale_tev,operator_strip_derivs,
    amplitude_survival_rank,unique_multi_derivative_projection,DERIVATIVE_KIND,
    round_trip_signature,model_strings,topology_key,_write_models,
    peak_memory_mib,
)
from neutrinomass.database.serialization import loads_completion
from neutrinomass.completions.fingerprints import (
    completion_fingerprint,species_model_fingerprint,propagator_model_fingerprint,
)


def artifact_path(directory,state,key="output"):
    return Path(directory)/state["artifacts"][key]["path"]


def summary_record(line):
    completion=loads_completion(line)
    fingerprint=completion_fingerprint(completion)
    projection=completion.lorentz_projection
    signature=(fingerprint,tuple(completion.momentum_contributions),
               tuple(projection) if projection is not None else None)
    return [json.dumps(dict(topology=topology_key(completion),
        topology_pair=[completion.topology,completion.canonical_topology],
        routed=bool(completion.momentum_contributions),models=model_strings(completion),
        species=repr(species_model_fingerprint(completion)),
        propagators=repr(propagator_model_fingerprint(completion)),
        fingerprint=repr(fingerprint),signature=repr(signature)),sort_keys=True)],{}


def digest_strings(values,*,sort=False):
    if sort:
        values=sorted(values)
    digest=hashlib.sha256()
    for i,value in enumerate(values):
        if i:
            digest.update(b"\n")
        digest.update(value.encode())
    return digest.hexdigest()


def run_operator(operator_name,historical_path,output_dir,contract,*,budget=None,scratch=None,
                 assignment_size=256,stage="all",unit=None,resume=False):
    """Execute stable work, then publish only a fully closed historical audit.

    ``unit`` is a zero-based index into the persisted generation plan. A
    selected generation unit includes decoded validation and local exact
    reduction; global merging always uses every unit in plan order.
    """
    started=time.monotonic()
    budget=budget or WorkBudget()
    root=Path(output_dir)
    historical_path=Path(historical_path)
    operator=operator_registry()[operator_name]
    order=operator_strip_derivs(operator.operator)["n_derivs"] if operator_kind(operator_name)==DERIVATIVE_KIND else 0
    contract=dict(contract,operator=operator_name,historical_sha256=digest_file(historical_path))
    plan=generation_plan(operator_name,assignment_size)
    plan_path=root/"plan.json"
    if plan_path.exists():
        if not resume:
            raise ValueError("Existing staged work requires --resume")
        if json.loads(plan_path.read_text()) != dict(contract=contract,units=plan):
            raise ValueError("Incompatible operator execution plan")
    else:
        atomic_json(plan_path,dict(contract=contract,units=plan))
    ranking=dict(representative_rank=amplitude_survival_rank,maximum_rank=1) if order and unique_multi_derivative_projection(operator) is not None else {}
    decoded_paths,local_paths,generation_states,decode_states=[],[],[],[]
    if unit is not None and (stage not in {"all","generation","decoded","local-reduction"} or unit not in range(len(plan))):
        raise ValueError("Unit selection only applies to a valid generation-plan unit")
    for i,selection in enumerate(plan):
        if unit is not None and unit != i:
            continue
        base=root/"units"/f"{i:08d}"
        generated=generate_unit(operator_name,selection,base/"generation",contract,budget=budget)
        generation_states.append(generated)
        if stage == "generation":
            continue
        decoded=map_records(output_paths(base/"generation",generated),base/"decoded",contract,
            decode_record,algorithm="decoded-validation-v2",budget=budget,scratch=scratch)
        decode_states.append((base/"decoded",decoded))
        decoded_paths.append(artifact_path(base/"decoded",decoded))
        if stage == "decoded":
            continue
        local=reduce_resumable([decoded_paths[-1]],base/"structural",contract,
            budget=budget,scratch=scratch,**ranking)
        local_paths.append(artifact_path(base/"structural",local))
    if stage in {"generation","decoded","local-reduction"} or unit is not None:
        return dict(operator=operator_name,stage=stage,unit=unit,complete=True,plan=str(plan_path))
    structural=reduce_resumable(local_paths,root/"structural",contract,budget=budget,scratch=scratch,**ranking)
    structural_path=artifact_path(root/"structural",structural)
    if stage == "structural":
        return structural
    amplitude=map_records([structural_path],root/"amplitudes",contract,amplitude_record,
        algorithm="amplitudes-v1",budget=budget,scratch=scratch)
    exact_path=artifact_path(root/"amplitudes",amplitude)
    if stage == "amplitudes":
        return amplitude
    historical=map_records([historical_path],root/"historical",contract,historical_transform(operator_name),
        algorithm="historical-normalization-v1",budget=budget,scratch=scratch)
    historical_raw=artifact_path(root/"historical",historical)
    historical_reduced=reduce_resumable([historical_raw],root/"historical-reduction",contract,budget=budget,scratch=scratch)
    historical_exact=artifact_path(root/"historical-reduction",historical_reduced)
    if stage == "historical":
        return historical_reduced
    coverage=historical_coverage([exact_path],[historical_exact],root/"coverage",contract,budget=budget,scratch=scratch)
    if coverage["progress"]["counts"].get("missing",0):
        raise RuntimeError("MISSING_HISTORICAL_CLASSES: production and publication blocked; see coverage evidence")
    if stage == "coverage":
        return coverage
    global_classes=reduce_resumable([exact_path],root/"global-lagrangians",contract,
        budget=budget,scratch=scratch,relation="global-lagrangian")
    # Compact summary records are restartable and avoid decoding every raw
    # completion again whenever final packaging is rerun.
    summaries={"generator":[artifact_path(directory,state,"evidence") for directory,state in decode_states]}
    for name,paths in {"exact_classes":[exact_path]}.items():
        result=map_records(paths,root/("summary-"+name),contract,summary_record,
            algorithm="summary-v1",budget=budget,scratch=scratch)
        summaries[name]=[artifact_path(root/("summary-"+name),result)]
    package_contract=dict(contract,algorithm="package-v1",inputs=input_contract([
        *decoded_paths,structural_path,exact_path,historical_exact,
        *(p for paths in summaries.values() for p in paths)]))
    with CheckpointStore(root/"publication",package_contract) as store:
        if store.state["complete"]:
            return json.loads(store.verify(store.state["artifacts"]["report"]).read_text())
        stats={}
        models=defaultdict(set)
        for name,paths in summaries.items():
            counts=Counter();topologies=Counter();species=set();propagators=set()
            signatures=hashlib.sha256()
            with tempfile.TemporaryDirectory(dir=scratch) as sort_directory:
                with sqlite3.connect(Path(sort_directory)/"fingerprints.sqlite3") as db:
                    db.execute("PRAGMA cache_size=-8192")
                    db.execute("PRAGMA temp_store=FILE")
                    db.execute("CREATE TABLE fingerprints (value TEXT)")
                    for line in lines_from_paths(paths):
                        data=json.loads(line)
                        if name == "generator" and data["code"] != "ACCEPTED":
                            continue
                        if counts["all"]:
                            signatures.update(b"\n")
                        counts["all"]+=1
                        counts["routed" if data["routed"] else "local"]+=1
                        topologies[data["topology"]]+=1
                        signatures.update(data["signature"].encode())
                        db.execute("INSERT INTO fingerprints VALUES (?)",(data["fingerprint"],))
                        if name == "exact_classes":
                            models[tuple(data["models"])].add(tuple(data["topology_pair"]))
                            species.add(data["species"]);propagators.add(data["propagators"])
                    digest=digest_strings(row[0] for row in db.execute("SELECT value FROM fingerprints ORDER BY value"))
            stats[name]=dict(counts={k:counts[k] for k in ("all","local","routed")},
                topologies=dict(sorted(topologies.items())),digest=digest,
                round_trip_digest=signatures.hexdigest(),species=len(species),propagators=len(propagators))
        invalid=[];historical_models=set()
        for line in lines_from_paths([artifact_path(root/"historical",historical,"evidence")]):
            data=json.loads(line);invalid.extend(data["invalid"])
            historical_models.add(tuple(data["representation_set"]))
        for line in lines_from_paths([artifact_path(root/"coverage",coverage)]):
            data=json.loads(line)
            if data["status"] == "invalid":
                invalid.append(data["evidence"])
        generation_rejections=Counter()
        for directory,state in decode_states:
            for line in lines_from_paths([artifact_path(directory,state,"evidence")]):
                data=json.loads(line)
                if data["code"] == "DECODED_VANISHING_VERTEX":
                    generation_rejections[data["topology"]]+=1
        rejected_amplitudes=Counter();sectors=Counter()
        for line in lines_from_paths([artifact_path(root/"amplitudes",amplitude,"evidence")]):
            data=json.loads(line);sectors[str(data["derivative_sectors"])]+=1
            if data["code"] == "AMPLITUDE_ZERO":
                rejected_amplitudes[topology_key(loads_completion(data["record"]))]+=1
        with tempfile.TemporaryDirectory(dir=scratch) as temporary:
            temporary=Path(temporary)
            # Concatenated gzip members preserve their original bytes and
            # remain an ordinary gzip JSONL stream for existing consumers.
            joined=temporary/"generator.jsonl.gz"
            with joined.open("wb") as dest:
                for path in decoded_paths:
                    with path.open("rb") as source:
                        shutil.copyfileobj(source,dest,1024**2)
            generator_artifact=store.publish_file(joined,suffix=".jsonl.gz")
            model_path=temporary/"models.jsonl"
            _write_models(model_path,models)
            model_artifact=store.publish_file(model_path,suffix=".jsonl")
        def reference(path):
            return dict(path=str(Path(path).resolve()),sha256=digest_file(path))
        physical=stats["exact_classes"];raw=stats["generator"]
        amplitude_report=dict(input_records=amplitude["progress"]["processed"],
            surviving_records=amplitude["progress"]["accepted"],rejected_records=amplitude["progress"]["rejected"],
            rejected_topologies=dict(sorted(rejected_amplitudes.items())),derivative_sectors=dict(sorted(sectors.items())),
            completion_digest=physical["digest"],source_sha256=digest_file(structural_path),destination_sha256=digest_file(exact_path))
        report=dict(schema_version=2,scale_unit="TeV",source_commit=contract.get("source_commit"),
            operator=operator_name,kind=operator_kind(operator_name),derivatives=order,
            hash_seed=contract.get("hash_seed",""),operator_scale_tev=operator_scale_tev(operator),
            wall_time_seconds=time.monotonic()-started,peak_memory_mib=peak_memory_mib(),
            stage_active_seconds={str(path.relative_to(root).parent):json.loads(path.read_text()).get("active_seconds",0)
                for path in sorted(root.rglob("checkpoint.json"))},
            records=dict(generator=raw["counts"]["all"],
                local=raw["counts"]["local"],routed=raw["counts"]["routed"]),
            provisional_records=sum(s["progress"]["records"] for s in generation_states),
            generation_rejections=dict(decoded_vanishing_uv_interactions=sum(generation_rejections.values()),topologies=dict(sorted(generation_rejections.items()))),
            structural_exact_classes=structural["progress"]["exact_classes"],exact_classes=physical["counts"],
            global_lagrangians=global_classes["progress"]["exact_classes"],democratic_models=len(models),
            species_models=physical["species"],propagator_models=physical["propagators"],
            historical=dict(records=historical["progress"]["processed"],
                normalized_records=historical["progress"]["accepted"],democratic_models=len(historical_models),
                classes=historical_reduced["progress"]["exact_classes"]+sum(len(json.loads(line)["invalid"]) for line in lines_from_paths([artifact_path(root/"historical",historical,"evidence")])),
                valid_classes=coverage["progress"]["counts"].get("reproduced",0),invalid_classes=len(invalid),invalid=invalid,
                amplitude_unsupported=[],reproduced=coverage["progress"]["counts"].get("reproduced",0),missing=0,
                exact_comparisons=coverage["progress"]["comparisons"],artifact=reference(historical_path)),
            topologies=dict(generator=raw["topologies"],exact_classes=physical["topologies"]),
            completion_digests=dict(generator=raw["digest"],structural_exact_classes=structural["progress"]["exact_completion_digest"],exact_classes=physical["digest"]),
            round_trip_digests=dict(generator=raw["round_trip_digest"],exact_classes=physical["round_trip_digest"]),
            amplitude_symmetrisation=amplitude_report,deduplication=structural["progress"],
            artifacts=dict(generator=reference(store.directory/generator_artifact["path"]),
                exact_classes=reference(exact_path),structural_exact_classes=reference(structural_path),
                democratic_models=reference(store.directory/model_artifact["path"]),
                historical_classes=reference(historical_exact),historical_coverage=reference(artifact_path(root/"coverage",coverage)),
                global_lagrangians=reference(artifact_path(root/"global-lagrangians",global_classes))),
            stage_manifest=dict(plan=reference(plan_path),generation_units=len(plan),
                topologies=sorted({u["topology"] for u in plan}),contract=contract))
        path=store.directory/".report.json"
        atomic_json(path,report)
        artifact=store.publish_file(path,suffix=".report.json")
        path.unlink()
        store.save(dict(published=True),artifacts=dict(report=artifact,generator=generator_artifact,models=model_artifact),complete=True)
        atomic_json(root/"report.json",report)
        return report
