from pathlib import Path
import pytest

from neutrinomass.database.pipeline import run_operator
from neutrinomass.database.checkpoints import WorkBudget,BudgetExpired
from rebuild_completion_database import stable_report_data


def test_complete_staged_operator_matches_recovered_run(tmp_path):
    historical=Path(__file__).resolve().parents[2]/"tests/fixtures/historical/op_1.dat"
    contract=dict(hash_seed="0",test="staged-pipeline")
    full=run_operator("1",historical,tmp_path/"full",contract,assignment_size=2)
    with pytest.raises(BudgetExpired):
        run_operator("1",historical,tmp_path/"resumed",contract,assignment_size=2,
                     budget=WorkBudget(fault_after=5))
    recovered=run_operator("1",historical,tmp_path/"resumed",contract,assignment_size=2,resume=True)
    assert stable_report_data(full) == stable_report_data(recovered)
    assert recovered["records"]["generator"] == 8
    assert recovered["exact_classes"]["all"] == recovered["global_lagrangians"] == 3
    assert recovered["historical"]["reproduced"] == 3
    assert recovered["historical"]["missing"] == 0
    assert run_operator("1",historical,tmp_path/"resumed",contract,assignment_size=2,resume=True) == recovered
