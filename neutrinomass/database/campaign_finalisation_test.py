from collections import Counter
import json
from pathlib import Path
import sqlite3
from types import SimpleNamespace

import pytest

from neutrinomass.completions.equivalence import field_label_parts
from neutrinomass.database.campaign_finalisation import (
    load_reports,verify_seeds,published_filter,build_species_index,species_key,
    known_tree_witness,filter_species_unit,SpeciesIndex,collect_species_results,global_reduction,model_comparison,
    historical_counting_replay,historical_interaction_key,publish_summary,
)
from neutrinomass.database.checkpoints import WorkBudget,BudgetExpired,digest_file
from neutrinomass.database.pipeline import run_operator,artifact_path
from neutrinomass.database.resumable_reduction import lines_from_paths
from neutrinomass.database.serialization import iter_completion_jsonl,dumps_completion
from neutrinomass.database.global_equivalence import species_descriptors


def test_complete_controls_have_exact_seed_proof_and_recoverable_filter_index(tmp_path):
    historical=Path(__file__).resolve().parents[2]/"tests/fixtures/historical/op_1.dat"
    entries={"1":{}}
    for seed in (0,1):
        root=tmp_path/f"seed-{seed}"
        run_operator("1",historical,root,dict(hash_seed=str(seed)))
        entries["1"][str(seed)]=dict(path=str(root/"report.json"),sha256=digest_file(root/"report.json"))
    with pytest.raises(ValueError,match="INVENTORY_INCOMPLETE"):
        load_reports(entries)
    reports=load_reports(entries,expected={"1"})
    assert verify_seeds(reports,tmp_path/"seed-agreement",{})[0]["physical_classes"]==3
    filtered=published_filter(reports,tmp_path/"published",{})
    assert filtered["progress"]["surviving_models"]==3
    with pytest.raises(BudgetExpired):
        build_species_index(reports,tmp_path/"species-index",{},budget=WorkBudget(fault_after=2))
    index=build_species_index(reports,tmp_path/"species-index",{})
    assert index["progress"]["species_multisets"]==3
    assert index["progress"]["physical_records"]==3
    path=artifact_path(tmp_path/"species-index",index,"database")
    reference=dict(path=str(path),sha256=digest_file(path))
    work=list(lines_from_paths([artifact_path(tmp_path/"species-index",index,"work")]))
    with SpeciesIndex(reference) as verified:
        for row in work:
            unit=json.loads(row)
            state=filter_species_unit(unit,verified,{"1":reports["1"]["0"]["operator_scale_tev"]},
                                      tmp_path/"species"/unit["identifier"],{})
            assert state["complete"] and state["progress"]["survives"]
    result=collect_species_results(tmp_path/"species-index",index,tmp_path/"species",tmp_path/"species-results",{},
        scales={"1":reports["1"]["0"]["operator_scale_tev"]})
    assert result["progress"]["accepted"]==3
    for mode,filter_path in (("all",None),("published",artifact_path(tmp_path/"published",filtered)),
                            ("species-aware",artifact_path(tmp_path/"species-results",result))):
        reduced=global_reduction(reports,tmp_path/"globals"/mode,{},mode=mode,filter_path=filter_path)
        assert reduced["progress"]["exact_classes"]==3
    old_models=tmp_path/"old-models.jsonl"
    rows=list(lines_from_paths([artifact_path(tmp_path/"published",filtered)]))
    old_models.write_text("\n".join(rows[:-1])+"\n")
    baseline=tmp_path/"baseline.json"
    baseline.write_text(json.dumps(dict(passed=True,per_operator_differences={},
        inputs={"1":dict(sha256=digest_file(historical))})))
    comparison=model_comparison(reports,artifact_path(tmp_path/"published",filtered),baseline,old_models,
                                tmp_path/"comparison",{})
    differences=list(map(json.loads,lines_from_paths([artifact_path(tmp_path/"comparison/differences",comparison)])))
    assert len(differences)==1
    assert differences[0]["mode"]=="published" and differences[0]["code"]=="MODEL_ADDED"
    with pytest.raises(BudgetExpired):
        historical_counting_replay(reports,artifact_path(tmp_path/"published",filtered),tmp_path/"historical-replay",{},
                                  budget=WorkBudget(fault_after=2))
    replay=historical_counting_replay(reports,artifact_path(tmp_path/"published",filtered),tmp_path/"historical-replay",{})
    assert replay["progress"]["interaction_multisets"]==3
    assert replay["progress"]["filtered_first_representative_lagrangians"]==3
    summary=publish_summary(entries,tmp_path,{},dict(path=str(baseline),sha256=digest_file(baseline)),expected={"1"})
    assert summary["complete_census"] is False
    assert summary["global_interaction_lagrangians"]==summary["global_species_multisets"]==3
    assert summary["filtered_species_aware"]["species_multisets"]==3
    assert summary["low_scale_genuineness_proved"] is False
    with pytest.raises(ValueError,match="INVENTORY_INCOMPLETE"):
        publish_summary(entries,tmp_path,{},dict(path=str(baseline),sha256=digest_file(baseline)))
    from scripts.finalise_erratum import run
    manifest=dict(entries=entries,producer_commit=None,
        baseline=dict(path=str(baseline),sha256=digest_file(baseline)),
        baseline_models=dict(path=str(old_models),sha256=digest_file(old_models)))
    with pytest.raises(BudgetExpired):
        run(manifest,tmp_path/"via-cli",{},expected={"1"},budget=WorkBudget(fault_after=4))
    resumed=run(manifest,tmp_path/"via-cli",{},expected={"1"})
    assert resumed["complete_census"] is False
    assert resumed["global_interaction_lagrangians"]==3
    assert resumed["filtered_species_aware"]==summary["filtered_species_aware"]
    # A completion certificate cannot mask damage to an underlying artifact.
    exact=Path(reports["1"]["0"]["artifacts"]["exact_classes"]["path"])
    exact.write_bytes(exact.read_bytes()[:-1])
    with pytest.raises(ValueError,match="HASH_MISMATCH"):
        load_reports(entries,expected={"1"})


def test_indexed_higher_tree_respects_copies_and_strictly_higher_scale():
    fixture=Path(__file__).resolve().parents[2]/"tests/fixtures/copy-mechanisms/0.jsonl"
    record=next(iter_completion_jsonl(fixture))
    descriptors=species_descriptors(record)
    repeated=next(value for value,count in Counter(descriptors.values()).items() if count==2)
    removed=next(label for label,description in descriptors.items() if description==repeated)
    one=SimpleNamespace(exotics={field for field in record.exotics
        if field_label_parts(field.label)[0]!=removed})
    with sqlite3.connect(":memory:") as db:
        db.execute("CREATE TABLE incidences(operator TEXT,key TEXT,scale REAL,payload TEXT)")
        db.execute("INSERT INTO incidences VALUES (?,?,?,?)",
                   (record.operator.name,species_key(record),100,dumps_completion(record)))
        assert known_tree_witness(record,10,db)["verified"]
        assert known_tree_witness(record,100,db) is None
        assert known_tree_witness(one,10,db) is None


def test_historical_interaction_replay_matches_literal_archive_head():
    import ast
    from neutrinomass.database.stages import historical_transform
    from neutrinomass.database.serialization import loads_completion
    from neutrinomass.utils.functions import conjugate_term
    fixture=Path(__file__).resolve().parents[2]/"tests/fixtures/historical/op_1.dat"
    for line in fixture.read_text().splitlines():
        if not line.strip():continue
        head=ast.literal_eval(line[len("LazyCompletion(head="):].split(", tail=",1)[0])
        expected=tuple(sorted(min(tuple(sorted(term)),conjugate_term(tuple(term))) for term in head["terms"]))
        records,_=historical_transform("1")(line)
        assert records
        assert all(historical_interaction_key(loads_completion(record))==expected for record in records)
