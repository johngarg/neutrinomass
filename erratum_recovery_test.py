from pathlib import Path
import pytest

from scripts.erratum_recovery import generation_replica,coverage_replica,reduction_replica,substitute
from neutrinomass.database.pipeline import run_operator
from neutrinomass.database.checkpoints import WorkBudget,BudgetExpired,digest_file


def test_recovered_generation_substitution_preserves_the_complete_scientific_result(tmp_path):
    historical=Path(__file__).parent/"tests/fixtures/historical"
    config=dict(historical=str(historical),assignment_size=2,recovery_operator="1")
    contract=dict(hash_seed="0",test="real-recovery-workflow")
    original=tmp_path/"original"
    run_operator("1",historical/"op_1.dat",original,contract,assignment_size=2)
    with pytest.raises(BudgetExpired):
        generation_replica(config,original,tmp_path/"generation",0,WorkBudget(fault_after=1),contract)
    recovered=generation_replica(config,original,tmp_path/"generation",0,WorkBudget(),contract)
    assert recovered["ordered_exact_agreement"]
    stage_contract=dict(contract,operator="1",historical_sha256=digest_file(historical/"op_1.dat"))
    with pytest.raises(BudgetExpired):
        coverage_replica(original,tmp_path/"coverage",0,WorkBudget(fault_after=1),stage_contract)
    covered=coverage_replica(original,tmp_path/"coverage",0,WorkBudget(),stage_contract)
    assert covered["exact_output_agreement"]
    with pytest.raises(BudgetExpired):
        reduction_replica(original,tmp_path/"reduction",WorkBudget(fault_after=1),stage_contract)
    reduced=reduction_replica(original,tmp_path/"reduction",WorkBudget(),stage_contract)
    assert reduced["ordered_exact_agreement"] and reduced["byte_identical"]
    assert reduced["records"]==3
    substituted=substitute(config,original,tmp_path/"substituted",tmp_path/"generation",0,WorkBudget(),contract)
    assert substituted["scientific_report_agreement"]
    assert substituted["ordered_exact_classes"]=={"exact_classes":3,"global_lagrangians":3}


def test_byte_identical_substitution_reuses_verified_stages_and_resumes_setup(tmp_path):
    historical=Path(__file__).parent/"tests/fixtures/historical"
    config=dict(historical=str(historical),assignment_size=2,recovery_operator="1")
    contract=dict(hash_seed="0",source_commit="tested",test="identical-input-substitution")
    original=tmp_path/"original"
    run_operator("1",historical/"op_1.dat",original,contract,assignment_size=2)
    before=(original/"report.json").read_bytes()
    with pytest.raises(BudgetExpired):
        substitute(config,original,tmp_path/"replica",original,0,WorkBudget(fault_after=2),contract)
    result=substitute(config,original,tmp_path/"replica",original,0,WorkBudget(),contract)
    assert result["method"]=="byte-identical-input-reuse"
    assert result["scientific_report_agreement"]
    assert result["ordered_exact_classes"]=={"exact_classes":3,"global_lagrangians":3}
    assert (original/"report.json").read_bytes()==before
    from scripts.validate_71p_gate import substitution_compatibility_gate,GateFailure
    substitution_compatibility_gate(result,"tested",operator="1")
    import json
    from neutrinomass.database.checkpoints import atomic_json,identity_hash
    checkpoint=Path(result["dependency_compatibility"]["path"])
    substituted=json.loads(Path(result["substituted_report"]["path"]).read_text())
    assert all(Path(artifact["path"]).is_relative_to(tmp_path/"replica") for artifact in substituted["artifacts"].values())
    assert substituted["recovery_substitution"]["retained_producer"]=="tested"
    proof=json.loads(checkpoint.read_text())
    proof["contract"]["consumer_inputs"].pop()
    proof["contract_sha256"]=identity_hash(proof["contract"])
    atomic_json(checkpoint,proof)
    result["dependency_compatibility"]["sha256"]=digest_file(checkpoint)
    with pytest.raises(GateFailure,match="CONSUMER_INPUTS_MISSING"):
        substitution_compatibility_gate(result,"tested",operator="1")


def test_substitution_rejects_overlapping_source_and_recovers_initial_lock(tmp_path):
    historical=Path(__file__).parent/"tests/fixtures/historical"
    config=dict(historical=str(historical),assignment_size=2,recovery_operator="1")
    contract=dict(hash_seed="0",test="initial-setup-crash")
    original=tmp_path/"original"
    run_operator("1",historical/"op_1.dat",original,contract,assignment_size=2)
    with pytest.raises(ValueError,match="OVERLAPS_SOURCE"):
        substitute(config,original,original,original,0,WorkBudget(),contract)
    replica=tmp_path/"replica"
    (replica/"substitution-setup").mkdir(parents=True)
    (replica/"substitution-setup/worker.lock").touch()
    result=substitute(config,original,replica,original,0,WorkBudget(),contract)
    assert result["method"]=="byte-identical-input-reuse"
