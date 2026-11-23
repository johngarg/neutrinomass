import json
from pathlib import Path
import pytest

from scripts.prepare_erratum_campaign import prepare,seal,reference
from neutrinomass.database.checkpoints import atomic_json
from neutrinomass.database.pipeline import run_operator


def test_execution_reuses_only_gate_verified_complete_stages_and_seals_both_seeds(tmp_path):
    historical=Path(__file__).parent/"tests/fixtures/historical"
    pilot=tmp_path/"pilot"
    reports=[]
    for seed in (0,1):
        root=tmp_path/"controls/tested"/f"seed-{seed}"/"1"
        run_operator("1",historical/"op_1.dat",root,dict(source_commit="tested",hash_seed=str(seed)),assignment_size=2)
        reports.append(reference(root/"report.json"))
    gate=tmp_path/"gate.json"
    atomic_json(gate,dict(passed=True,producer_commit="tested",reports=[],controls=reports))
    ledger=tmp_path/"ledger.json";atomic_json(ledger,dict(production_gate=reference(gate)))
    config=tmp_path/"pilot-config.json"
    atomic_json(config,dict(source_commit="tested",root=str(pilot),ledger=str(ledger),historical=str(historical)))
    baseline=tmp_path/"baseline.json"
    atomic_json(baseline,dict(passed=True,inputs={"1":reference(historical/"op_1.dat")}))
    models=tmp_path/"models.jsonl";models.write_text("baseline\n")
    args=(config,reference(baseline),reference(models),tmp_path/"production",dict(commit="validator"))
    execution=prepare(*args,expected={"1"},reused={"1"})
    assert prepare(*args,expected={"1"},reused={"1"})==execution
    result=json.loads(execution.read_text())
    assert result["production_config"]["operators"]==[]
    assert result["producer_commit"]=="tested"
    reuse=json.loads((execution.parent/"reuse-manifest.json").read_text())
    assert reuse["operators"]["1"]["0"]["stages"]
    assert all(row["reuse"] for row in reuse["operators"]["1"]["0"]["stages"])
    sealed=seal(execution,tmp_path/"finalisation-inputs.json")
    assert set(sealed["entries"]["1"])=={"0","1"}
    assert sealed["entries"]["1"]["0"]==reports[0]
    # A formerly valid report alone cannot establish reusable stage data.
    checkpoint=tmp_path/"controls/tested/seed-0/1/units/00000000/generation/checkpoint.json"
    checkpoint.unlink()
    with pytest.raises(ValueError,match="REQUIRED_UNIT_MISSING"):
        prepare(*args,expected={"1"},reused={"1"})


def test_production_execution_plan_requires_intact_gate_and_exact_inventory(tmp_path):
    config=tmp_path/"config.json";ledger=tmp_path/"ledger.json"
    atomic_json(config,dict(source_commit="tested",ledger=str(ledger)))
    atomic_json(ledger,dict(production_gate=None))
    with pytest.raises(ValueError,match="PASSED_71P_GATE"):
        prepare(config,{}, {},tmp_path/"out",{})
    gate=tmp_path/"gate.json";atomic_json(gate,dict(passed=True,producer_commit="tested"))
    atomic_json(ledger,dict(production_gate=reference(gate)))
    gate.write_text("corrupt")
    with pytest.raises(ValueError,match="REFERENCE_CHANGED"):
        prepare(config,{}, {},tmp_path/"out",{})


def test_sealing_requires_completed_queue_publication_not_merely_a_report(tmp_path):
    gate=tmp_path/"gate.json";atomic_json(gate,dict(passed=True,producer_commit="tested"))
    ledger=tmp_path/"ledger.json";atomic_json(ledger,dict(production_gate=reference(gate)))
    report=tmp_path/"report.json";atomic_json(report,{})
    execution=tmp_path/"execution.json"
    atomic_json(execution,dict(production_config=dict(ledger=str(ledger),source_commit="tested",queue=str(tmp_path/"queue.sqlite3")),
        entries={"1":{"0":dict(action="recompute",report_path=str(report))}}))
    with pytest.raises(ValueError,match="OWNING_QUEUE"):
        seal(execution,tmp_path/"final.json")
    assert not (tmp_path/"final.json").exists()
