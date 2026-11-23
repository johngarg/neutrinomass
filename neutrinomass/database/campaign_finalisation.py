"""Independently runnable finalisation of a complete two-seed census.

Publication requires the complete operator inventory. Exact seed comparison,
cross-operator reduction and candidate mechanism searches have their own
durable checkpoints; no scheduler dependency chain is needed to restart them.
"""

from collections import defaultdict
import ast
from itertools import combinations
import json
from pathlib import Path
import sqlite3

from neutrinomass.database.checkpoints import (
    CheckpointStore,WorkBudget,atomic_json,canonical_json,digest_file,identity_hash,
)
from neutrinomass.database.filtering import (
    filter_published_registry,published_model_key,multiset_contains,ONE_LOOP_SCALE_TEV,
)
from neutrinomass.database.global_equivalence import species_descriptors
from neutrinomass.database.mechanism_checkpoints import tree_search,loop_search,read_witness
from neutrinomass.database.pipeline import artifact_path
from neutrinomass.database.rebuild import operator_registry,model_strings
from neutrinomass.database.resumable_reduction import reduce_resumable,lines_from_paths,input_contract
from neutrinomass.database.serialization import loads_completion
from neutrinomass.database.stages import map_records


def verify_reference(reference):
    path=Path(reference["path"])
    if not path.is_file() or digest_file(path)!=reference["sha256"]:
        raise ValueError("CAMPAIGN_ARTIFACT_HASH_MISMATCH: "+str(path))
    return path


def load_reports(entries,*,expected=None,verify_artifacts=True):
    """Validate coverage and all referenced bytes before accepting any totals."""
    expected=set(operator_registry()) if expected is None else set(expected)
    if set(entries)!=expected:
        raise ValueError("CAMPAIGN_OPERATOR_INVENTORY_INCOMPLETE")
    reports={}
    for operator in sorted(entries):
        if set(entries[operator])!={"0","1"}:
            raise ValueError("CAMPAIGN_SEED_INVENTORY_INCOMPLETE")
        reports[operator]={}
        for seed,reference in entries[operator].items():
            report=json.loads(verify_reference(reference).read_text())
            if (report.get("schema_version")!=2 or report.get("scale_unit")!="TeV"
                    or report["operator"]!=operator or str(report["hash_seed"])!=seed):
                raise ValueError("CAMPAIGN_REPORT_IDENTITY_MISMATCH")
            if report["historical"]["missing"] or report["historical"]["amplitude_unsupported"]:
                raise ValueError("CAMPAIGN_UNRESOLVED_SCIENTIFIC_DIAGNOSTIC")
            if verify_artifacts:
                for artifact in report["artifacts"].values():
                    verify_reference(artifact)
                verify_reference(report["historical"]["artifact"])
            plan=json.loads(verify_reference(report["stage_manifest"]["plan"]).read_text())
            if len(plan["units"])!=report["stage_manifest"]["generation_units"]:
                raise ValueError("CAMPAIGN_GENERATION_PLAN_INCOMPLETE")
            reports[operator][seed]=report
    return reports


def verify_seeds(reports,directory,contract,*,budget=None,scratch=None):
    from rebuild_completion_database import stable_report_data
    budget=budget or WorkBudget()
    evidence=[]
    for operator,seeds in sorted(reports.items()):
        a,b=seeds["0"],seeds["1"]
        if stable_report_data(a)!=stable_report_data(b):
            raise ValueError("CAMPAIGN_SCIENTIFIC_SEED_DISAGREEMENT: "+operator)
        union=reduce_resumable([r["artifacts"]["exact_classes"]["path"] for r in (a,b)],
            Path(directory)/operator,contract,budget=budget,scratch=scratch)
        if union["progress"]["exact_classes"]!=a["exact_classes"]["all"]:
            raise ValueError("CAMPAIGN_EXACT_SEED_CLASSES_DIFFER: "+operator)
        models=[{published_model_key(json.loads(line)["fields"])
                 for line in lines_from_paths([r["artifacts"]["democratic_models"]["path"]])}
                for r in (a,b)]
        if models[0]!=models[1]:
            raise ValueError("CAMPAIGN_SEED_MODEL_SETS_DIFFER: "+operator)
        evidence.append(dict(operator=operator,physical_classes=union["progress"]["exact_classes"],
            checkpoint=dict(path=str(Path(directory)/operator/"checkpoint.json"),
                            sha256=digest_file(Path(directory)/operator/"checkpoint.json"))))
    return evidence


def published_filter(reports,directory,contract):
    registry={};scales={};paths=[]
    for operator,seeds in sorted(reports.items()):
        report=seeds["0"];path=report["artifacts"]["democratic_models"]["path"]
        paths.append(path)
        registry[operator]={tuple(json.loads(line)["fields"]):None for line in lines_from_paths([path])}
        scales[operator]=report["operator_scale_tev"]
    identity=dict(contract,algorithm="campaign-published-filter-v1",inputs=input_contract(paths),scales=scales)
    with CheckpointStore(directory,identity) as store:
        if store.state["complete"]:
            return store.state
        result=filter_published_registry(registry,scales)
        survivors=store.publish_lines(canonical_json(dict(operator=op,fields=fields))
            for op,models in sorted(result["survivors"].items()) for fields in sorted(models))
        witnesses=store.publish_lines(canonical_json(row) for row in result["witnesses"])
        per_operator={op:dict(input=len(registry[op]),survivors=len(result["survivors"][op]),
            removed_mass=result["removed_by_mass"].get(op,0),
            removed_loop=result["removed_by_one_loop_weinberg"].get(op,0)) for op in registry}
        store.save(dict(mode="published",scale_unit="TeV",operators=per_operator,
            input_models=result["distinct_input_models"],surviving_models=result["distinct_surviving_models"]),
            artifacts=dict(output=survivors,witnesses=witnesses),complete=True)
        return store.state


def species_key(completion):
    return canonical_json(sorted(species_descriptors(completion).values(),key=repr))


def build_species_index(reports,directory,contract,*,budget=None,scratch=None):
    """Store one concrete representative per actual species multiset on disk."""
    budget=budget or WorkBudget()
    operators=sorted(reports)
    paths=[reports[op]["0"]["artifacts"]["exact_classes"]["path"] for op in operators]
    scales={op:reports[op]["0"]["operator_scale_tev"] for op in operators}
    identity=dict(contract,algorithm="campaign-species-index-v1",operators=operators,
                  inputs=input_contract(paths),scales=scales)
    with CheckpointStore(directory,identity) as store:
        if store.state["complete"]:
            return store.state
        with store.database(scratch) as db:
            db.executescript("""
                CREATE TABLE IF NOT EXISTS species (key TEXT PRIMARY KEY,payload TEXT NOT NULL);
                CREATE TABLE IF NOT EXISTS incidences (operator TEXT,key TEXT,scale REAL,payload TEXT,
                    PRIMARY KEY(operator,key));
                CREATE INDEX IF NOT EXISTS by_species ON incidences(key,scale);
            """)
            progress=dict(store.state["progress"] or dict(operator_index=0,record=0,physical_records=0))
            def save():store.snapshot(db,progress)
            for i in range(progress["operator_index"],len(operators)):
                op=operators[i]
                for ordinal,line in enumerate(lines_from_paths([paths[i]]),1):
                    if ordinal<=progress["record"]:continue
                    key=species_key(loads_completion(line))
                    db.execute("INSERT OR IGNORE INTO species VALUES (?,?)",(key,line))
                    db.execute("INSERT OR IGNORE INTO incidences VALUES (?,?,?,?)",(op,key,scales[op],line))
                    progress.update(record=ordinal,physical_records=progress["physical_records"]+1)
                    budget.boundary(save,progress=dict(stage="species-index",operator=op,**progress))
                progress.update(operator_index=i+1,record=0)
            progress["species_multisets"]=db.execute("SELECT COUNT(*) FROM species").fetchone()[0]
            progress["operator_species_incidences"]=db.execute("SELECT COUNT(*) FROM incidences").fetchone()[0]
            work=store.publish_lines(canonical_json(dict(species=key,scale_tev=scale,
                identifier=identity_hash(dict(species=key,scale_tev=scale))))
                for key,scale in db.execute("SELECT DISTINCT key,scale FROM incidences ORDER BY key,scale DESC"))
            store.save(progress,artifacts={"work":work})
            store.snapshot(db,progress,complete=True)
            return store.state


def _frozen(value):
    return tuple(_frozen(x) for x in value) if isinstance(value,list) else value


def known_tree_witness(candidate,scale,index):
    """Use a proved census mechanism only when the candidate admits its copies."""
    from neutrinomass.completions.amplitudes import audit_amplitude_symmetrisation
    from neutrinomass.database.rebuild import validate_completion
    available=tuple(sorted(species_descriptors(candidate).values(),key=repr))
    for size in range(1,len(available)+1):
        for subset in dict.fromkeys(combinations(available,size)):
            key=canonical_json(subset)
            found=index.execute("""SELECT i.operator,i.scale,i.payload FROM incidences i
                WHERE i.key=? AND i.scale>?
                ORDER BY i.scale DESC,i.operator LIMIT 1""",(key,scale)).fetchone()
            if found is None:continue
            operator,mechanism_scale,payload=found
            record=loads_completion(payload)
            if record.operator.name!=operator:
                raise ValueError("CENSUS_MECHANISM_OPERATOR_MISMATCH")
            validate_completion(record)
            audit=audit_amplitude_symmetrisation(record)
            if audit.status!="nonzero":
                raise ValueError("CENSUS_MECHANISM_WITNESS_NOT_PHYSICAL")
            required=tuple(sorted(species_descriptors(record).values(),key=repr))
            if not multiset_contains(available,required):
                raise ValueError("CENSUS_MECHANISM_REQUIRES_ABSENT_COPY")
            return dict(code="HIGHER_SCALE_CENSUS_TREE",verified=True,
                operator=record.operator.name,scale_tev=mechanism_scale,required_species=required,
                species_assignment=species_descriptors(record),record=payload,
                amplitude_witness=audit.witness)
    return None


class SpeciesIndex:
    """Verify the immutable snapshot once per worker, then open it read-only."""
    def __init__(self,reference):
        self.reference=reference
        self.connection=None

    def __enter__(self):
        path=verify_reference(self.reference)
        self.connection=sqlite3.connect(path.resolve().as_uri()+"?mode=ro",uri=True)
        return self

    def __exit__(self,*args):
        self.connection.close()
        self.connection=None


def filter_species_unit(work,index,scales,directory,contract,*,budget=None):
    budget=budget or WorkBudget()
    if index.connection is None:
        raise ValueError("FILTER_REQUIRES_VERIFIED_SPECIES_INDEX")
    identity=dict(contract,algorithm="species-aware-filter-unit-v1",work=work,
        index_sha256=index.reference["sha256"],scales=scales)
    directory=Path(directory)
    with CheckpointStore(directory,identity) as store:
        if store.state["complete"]:return store.state
        row=index.connection.execute("SELECT payload FROM species WHERE key=?",(work["species"],)).fetchone()
        if row is None:raise ValueError("FILTER_CANDIDATE_NOT_IN_VERIFIED_INDEX")
        candidate=loads_completion(row[0])
        witness=known_tree_witness(candidate,work["scale_tev"],index.connection)
        if witness is None:
            higher=tuple(sorted((op for op in scales if scales[op]>work["scale_tev"]),
                                key=lambda op:(-scales[op],op)))
            state=tree_search(candidate,higher,directory/"tree",contract,budget=budget)
            witness=read_witness(directory/"tree",state)
        if witness is None and work["scale_tev"]<ONE_LOOP_SCALE_TEV:
            state=loop_search(candidate,directory/"loops",contract,budget=budget)
            witness=read_witness(directory/"loops",state)
        if witness is not None and not multiset_contains(
                _frozen(json.loads(work["species"])),_frozen(witness["required_species"])):
            raise ValueError("FILTER_WITNESS_REQUIRES_ABSENT_COPY")
        result=dict(work,mode="species-aware",scale_unit="TeV",survives=witness is None,
            witness=witness,mechanism_domain="published scalar/fermion trees and heavy T1/T3 Weinberg loops",
            low_scale_genuineness_proved=False)
        output=store.publish_lines([canonical_json(result)])
        store.save(dict(survives=witness is None),artifacts={"output":output},complete=True)
        return store.state


def collect_species_results(index_directory,index_state,unit_root,directory,contract,
                            *,scales,budget=None,scratch=None):
    work_path=artifact_path(index_directory,index_state,"work")
    index_hash=index_state["artifacts"]["database"]["sha256"]
    def collect(line):
        work=json.loads(line)
        identity=dict(contract,algorithm="species-aware-filter-unit-v1",work=work,
                      index_sha256=index_hash,scales=scales)
        with CheckpointStore(Path(unit_root)/work["identifier"],identity) as unit:
            if not unit.state["complete"]:
                raise ValueError("SPECIES_FILTER_WORK_INCOMPLETE: "+work["identifier"])
            rows=list(lines_from_paths([unit.verify(unit.state["artifacts"]["output"])]))
            if len(rows)!=1 or json.loads(rows[0])["identifier"]!=work["identifier"]:
                raise ValueError("SPECIES_FILTER_RESULT_IDENTITY_MISMATCH")
            return rows,{}
    return map_records([work_path],directory,dict(contract,index_sha256=index_hash,scales=scales),collect,
                       algorithm="species-result-collection-v1",budget=budget,scratch=scratch)


def global_reduction(reports,directory,contract,*,mode="all",filter_path=None,budget=None,scratch=None):
    """Union across operators using exact interaction-multiset equivalence."""
    if mode not in {"all","published","species-aware"}:
        raise ValueError("UNKNOWN_GLOBAL_FILTER_MODE")
    budget=budget or WorkBudget()
    directory=Path(directory)
    survivors=defaultdict(set)
    admitted=set()
    identity=dict(contract,filter_mode=mode)
    if mode!="all":
        if filter_path is None:raise ValueError("GLOBAL_FILTER_INPUT_REQUIRED")
        identity["filter_input"]=input_contract([filter_path])
        for line in lines_from_paths([filter_path]):
            row=json.loads(line)
            if mode=="published":survivors[row["operator"]].add(published_model_key(row["fields"]))
            elif row["survives"]:admitted.add((row["species"],row["scale_tev"]))
    paths=[]
    for operator,seeds in sorted(reports.items()):
        report=seeds["0"];source=report["artifacts"]["global_lagrangians"]["path"]
        if mode=="all":
            paths.append(source)
            continue
        def select(line):
            record=loads_completion(line)
            keep=(published_model_key(model_strings(record)) in survivors[operator] if mode=="published"
                  else (species_key(record),report["operator_scale_tev"]) in admitted)
            return ([line] if keep else []),{}
        selection=map_records([source],directory/"selected"/operator,identity,select,
            algorithm="global-filter-selection-v1",budget=budget,scratch=scratch)
        paths.append(artifact_path(directory/"selected"/operator,selection))
    return reduce_resumable(paths,directory/"reduction",identity,
        relation="global-lagrangian",budget=budget,scratch=scratch)


def model_comparison(reports,published_path,baseline_path,baseline_models,directory,contract,
                     *,budget=None,scratch=None):
    baseline=json.loads(Path(baseline_path).read_text())
    if not baseline["passed"] or baseline["per_operator_differences"]:
        raise ValueError("HISTORICAL_BASELINE_NOT_REPRODUCED")
    old_filtered=defaultdict(dict);new_filtered=defaultdict(dict)
    for path,registry in ((baseline_models,old_filtered),(published_path,new_filtered)):
        for line in lines_from_paths([path]):
            row=json.loads(line)
            registry[row["operator"]][published_model_key(row["fields"])]=tuple(row["fields"])
    paths=[baseline_path,baseline_models,published_path]
    for operator,seeds in sorted(reports.items()):
        report=seeds["0"]
        if baseline["inputs"][operator]["sha256"]!=report["historical"]["artifact"]["sha256"]:
            raise ValueError("HISTORICAL_BASELINE_INPUT_MISMATCH")
        paths.extend([report["historical"]["artifact"]["path"],report["artifacts"]["democratic_models"]["path"]])
    identity=dict(contract,algorithm="campaign-model-comparison-v1",inputs=input_contract(paths))
    directory=Path(directory)
    with CheckpointStore(directory/"inputs",identity) as store:
        if not store.state["complete"]:
            output=store.publish_lines(canonical_json(dict(operator=op)) for op in sorted(reports))
            store.save(dict(operators=len(reports)),artifacts={"output":output},complete=True)
        source=artifact_path(directory/"inputs",store.state)
    def compare(line):
        operator=json.loads(line)["operator"];report=reports[operator]["0"]
        old={}
        for historical in lines_from_paths([report["historical"]["artifact"]["path"]]):
            prefix="LazyCompletion(head="
            if not historical.startswith(prefix):raise ValueError("UNKNOWN_HISTORICAL_HEADER")
            head=ast.literal_eval(historical[len(prefix):].split(", tail=",1)[0])
            if head["operator_name"]!=operator:raise ValueError("HISTORICAL_OPERATOR_MISMATCH")
            fields=tuple(sorted(set(head["quantum_numbers"])))
            old[published_model_key(fields)]=fields
        new={published_model_key(row["fields"]):tuple(row["fields"])
             for row in map(json.loads,lines_from_paths([report["artifacts"]["democratic_models"]["path"]]))}
        output=[];counts={}
        for mode,before,after in (("unfiltered",old,new),("published",old_filtered[operator],new_filtered[operator])):
            added=set(after)-set(before);removed=set(before)-set(after)
            counts[mode]=dict(historical=len(before),current=len(after),added=len(added),removed=len(removed))
            for code,keys,lookup in (("MODEL_ADDED",added,after),("MODEL_REMOVED",removed,before)):
                output.extend(canonical_json(dict(operator=operator,mode=mode,code=code,key=key,fields=lookup[key]))
                              for key in sorted(keys))
        return output,dict(operator=operator,counts=counts)
    return map_records([source],directory/"differences",identity,compare,
        algorithm="operator-model-differences-v1",budget=budget,scratch=scratch)


def historical_interaction_key(completion):
    """Replay the published contraction-blind key, separately from exact equality."""
    from neutrinomass.utils.functions import stringify_qns,conjugate_term
    terms=[tuple(sorted(stringify_qns(field) for field in term.fields)) for term in completion.terms]
    return tuple(sorted(min(term,conjugate_term(term)) for term in terms))


def historical_counting_replay(reports,published_path,directory,contract,*,budget=None,scratch=None):
    budget=budget or WorkBudget()
    surviving=defaultdict(set)
    for line in lines_from_paths([published_path]):
        row=json.loads(line)
        surviving[row["operator"]].add(published_model_key(row["fields"]))
    operators=sorted(reports)
    paths=[reports[op]["0"]["artifacts"]["exact_classes"]["path"] for op in operators]
    identity=dict(contract,algorithm="historical-counting-replay-v1",operators=operators,
                  inputs=input_contract(paths),filter_input=input_contract([published_path]))
    with CheckpointStore(directory,identity) as store:
        if store.state["complete"]:return store.state
        with store.database(scratch) as db:
            db.executescript("""
                CREATE TABLE IF NOT EXISTS terms (key TEXT PRIMARY KEY);
                CREATE TABLE IF NOT EXISTS first_models (operator TEXT,model TEXT,terms TEXT,survives INTEGER,
                    PRIMARY KEY(operator,model));
            """)
            progress=dict(store.state["progress"] or dict(operator_index=0,record=0,physical_records=0))
            def save():store.snapshot(db,progress)
            for i in range(progress["operator_index"],len(operators)):
                op=operators[i]
                for ordinal,line in enumerate(lines_from_paths([paths[i]]),1):
                    if ordinal<=progress["record"]:continue
                    record=loads_completion(line)
                    key=canonical_json(historical_interaction_key(record))
                    fields=model_strings(record)
                    db.execute("INSERT OR IGNORE INTO terms VALUES (?)",(key,))
                    db.execute("INSERT OR IGNORE INTO first_models VALUES (?,?,?,?)",
                        (op,canonical_json(fields),key,int(published_model_key(fields) in surviving[op])))
                    progress.update(record=ordinal,physical_records=progress["physical_records"]+1)
                    budget.boundary(save,progress=dict(stage="historical-counting-replay",operator=op,**progress))
                progress.update(operator_index=i+1,record=0)
            progress.update(interaction_multisets=db.execute("SELECT COUNT(*) FROM terms").fetchone()[0],
                filtered_first_representative_lagrangians=db.execute(
                    "SELECT COUNT(DISTINCT terms) FROM first_models WHERE survives=1").fetchone()[0])
            store.snapshot(db,progress,complete=True)
            return store.state


def completed_stage(directory,contract=None):
    path=Path(directory)/"checkpoint.json"
    if not path.is_file():raise ValueError("FINALISATION_STAGE_MISSING: "+str(directory))
    original=json.loads(path.read_text())
    if contract is not None and any(original["contract"].get(key)!=value for key,value in contract.items()):
        raise ValueError("FINALISATION_STAGE_PRODUCER_MISMATCH")
    with CheckpointStore(directory,original["contract"]) as checked:
        if not checked.state["complete"]:
            raise ValueError("FINALISATION_STAGE_INCOMPLETE: "+str(directory))
        return checked.state


def publish_summary(entries,directory,contract,baseline_reference,*,expected=None,budget=None,scratch=None):
    """Publish complete totals only after every prerequisite is verified."""
    directory=Path(directory)
    reports=load_reports(entries,expected=expected)
    full=set(reports)==set(operator_registry()) and len(reports)==243
    if not full and expected is None:
        raise ValueError("COMPLETE_CENSUS_REQUIRES_243_OPERATORS")
    agreement=verify_seeds(reports,directory/"seed-agreement",contract,budget=budget,scratch=scratch)
    published=completed_stage(directory/"published",contract)
    index=completed_stage(directory/"species-index",contract)
    species=completed_stage(directory/"species-results",contract)
    comparison=completed_stage(directory/"comparison/differences",contract)
    replay=completed_stage(directory/"historical-replay",contract)
    global_counts={mode:completed_stage(directory/"globals"/mode/"reduction",contract)
                   for mode in ("all","published","species-aware")}
    baseline=json.loads(verify_reference(baseline_reference).read_text())
    if (not baseline["passed"] or baseline["per_operator_differences"] or
            (full and (baseline["literal_unfiltered_models"],baseline["filtered_models"],
             baseline["historical_interaction_multisets"],baseline["historical_filtered_first_representative_lagrangians"])
             !=(141989,11216,430810,11483))):
        raise ValueError("FINALISATION_BASELINE_INVALID")
    results={}
    for line in lines_from_paths([artifact_path(directory/"species-results",species)]):
        row=json.loads(line)
        key=(row["species"],row["scale_tev"])
        if key in results:raise ValueError("DUPLICATED_SPECIES_FILTER_RESULT")
        results[key]=row["survives"]
    species_counts=defaultdict(lambda:dict(input=0,survivors=0))
    index_reference=dict(path=str(artifact_path(directory/"species-index",index,"database")),
                         sha256=index["artifacts"]["database"]["sha256"])
    expected_work=set()
    with SpeciesIndex(index_reference) as verified:
        for operator,key,scale in verified.connection.execute("SELECT operator,key,scale FROM incidences"):
            expected_work.add((key,scale))
            if (key,scale) not in results:raise ValueError("SPECIES_FILTER_WORK_INCOMPLETE")
            species_counts[operator]["input"]+=1
            species_counts[operator]["survivors"]+=results[key,scale]
    if set(results)!=expected_work:raise ValueError("SPECIES_FILTER_WORK_INVENTORY_MISMATCH")
    model_changes={json.loads(line)["operator"]:json.loads(line)["counts"] for line in
        lines_from_paths([artifact_path(directory/"comparison/differences",comparison,"evidence")])}
    rows={}
    for operator,seeds in sorted(reports.items()):
        report=seeds["0"]
        rows[operator]=dict(generated_records=report["records"]["generator"],
            structural_classes=report["structural_exact_classes"],physical_completion_classes=report["exact_classes"]["all"],
            interaction_lagrangians=report["global_lagrangians"],representation_sets=report["democratic_models"],
            species_multisets=species_counts[operator]["input"],species_aware_survivors=species_counts[operator]["survivors"],
            published_filter=published["progress"]["operators"][operator],model_changes=model_changes[operator],
            historical_reproduced=report["historical"]["reproduced"],
            historical_invalid=report["historical"]["invalid_classes"],historical_missing=report["historical"]["missing"])
    stages={"published":directory/"published","species_index":directory/"species-index",
            "species_results":directory/"species-results","comparison":directory/"comparison/differences",
            "historical_replay":directory/"historical-replay"}
    stages.update({"global_"+mode:directory/"globals"/mode/"reduction" for mode in global_counts})
    dependencies={key:dict(path=str(path/"checkpoint.json"),sha256=digest_file(path/"checkpoint.json"))
                  for key,path in stages.items()}
    identity=dict(contract,algorithm="complete-census-publication-v1",reports=entries,
                  dependencies=dependencies,baseline=baseline_reference)
    with CheckpointStore(directory/"publication",identity) as store:
        if store.state["complete"]:
            return json.loads(store.verify(store.state["artifacts"]["report"]).read_text())
        result=dict(schema_version=2,scale_unit="TeV",complete_census=full,operator_count=len(reports),seeds=[0,1],
            scientific_versions=contract,source_reports=entries,exact_seed_agreement=agreement,
            generated_records=sum(row["generated_records"] for row in rows.values()),
            physical_completion_operator_incidences=sum(row["physical_completion_classes"] for row in rows.values()),
            global_interaction_lagrangians=global_counts["all"]["progress"]["exact_classes"],
            global_representation_sets=published["progress"]["input_models"],
            global_species_multisets=index["progress"]["species_multisets"],
            filtered_published=dict(representation_sets=published["progress"]["surviving_models"],
                global_interaction_lagrangians=global_counts["published"]["progress"]["exact_classes"]),
            filtered_species_aware=dict(species_multisets=len({key for (key,scale),survives in results.items() if survives}),
                global_interaction_lagrangians=global_counts["species-aware"]["progress"]["exact_classes"]),
            historical_counting_replay=dict(interaction_multisets=replay["progress"]["interaction_multisets"],
                filtered_first_representative_lagrangians=replay["progress"]["filtered_first_representative_lagrangians"]),
            historical_baseline=baseline_reference,operators=rows,dependencies=dependencies,
            assumptions=["generic allowed couplings","candidate baryon assignments","no additional imposed symmetries"],
            mechanism_domain="published scalar/fermion trees and heavy T1/T3 Weinberg loops",
            low_scale_genuineness_proved=False,
            limitations=["The vector extension is excluded.",
                         "Passing these finite mechanism filters does not prove low-scale genuineness."])
        temporary=store.directory/".report.json"
        atomic_json(temporary,result)
        artifact=store.publish_file(temporary,suffix=".report.json")
        temporary.unlink()
        store.save(dict(complete_census=full),artifacts={"report":artifact},complete=True)
        atomic_json(directory/"report.json",result)
        return result
