from pathlib import Path
import pytest

from scripts.erratum_recovery import generation_replica,coverage_replica,substitute
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
    substituted=substitute(config,original,tmp_path/"substituted",tmp_path/"generation",0,WorkBudget(),contract)
    assert substituted["scientific_report_agreement"]
    assert substituted["ordered_exact_classes"]=={"exact_classes":3,"global_lagrangians":3}
