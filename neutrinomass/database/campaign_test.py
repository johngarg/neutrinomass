import pytest
import json
from neutrinomass.database.campaign import CampaignLedger, CampaignBudgetError
from neutrinomass.database.checkpoints import digest_file


def test_reservations_include_failed_allocations_and_inflight_maximum(tmp_path):
    ledger=CampaignLedger(tmp_path/"ledger.json")
    for i in range(4):
        ledger.reserve(str(i))
    with pytest.raises(CampaignBudgetError,match="concurrency"):
        ledger.reserve("fifth")
    ledger.finished("0",elapsed_seconds=1800,allocated_cpus=1,status="FAILED",failure_class="infrastructure")
    ledger.reserve("retry",retry_of="0")
    assert ledger.summary()["charged_core_hours"] == 8.5
    with pytest.raises(CampaignBudgetError,match="gate"):
        ledger.reserve("production",phase="production")


def test_initial_pilot_budget_is_not_silently_extended(tmp_path):
    ledger=CampaignLedger(tmp_path/"ledger.json")
    for i in range(100):
        ledger.reserve(str(i))
        ledger.finished(str(i),elapsed_seconds=7200,allocated_cpus=1,status="COMPLETED")
    with pytest.raises(CampaignBudgetError,match="200"):
        ledger.reserve("over")
    assert ledger.summary()["remaining_core_hours"] == 2800


def test_scientific_failure_and_retry_limit_require_diagnosis(tmp_path):
    ledger=CampaignLedger(tmp_path/"ledger.json")
    ledger.reserve("science")
    ledger.finished("science",elapsed_seconds=3,allocated_cpus=1,status="FAILED",failure_class="scientific")
    with pytest.raises(ValueError,match="diagnosed"):
        ledger.reserve("bad-retry",retry_of="science")
    for i in range(3):
        ledger.reserve(str(i),retry_of=str(i-1) if i else None)
        ledger.finished(str(i),elapsed_seconds=3,allocated_cpus=1,status="FAILED",failure_class="infrastructure")
    with pytest.raises(ValueError,match="exhausted"):
        ledger.reserve("third-retry",retry_of="2")


def test_production_revalidates_gate_bytes_before_every_reservation(tmp_path):
    ledger=CampaignLedger(tmp_path/"ledger.json")
    gate=tmp_path/"gate.json"
    gate.write_text(json.dumps(dict(passed=True,producer_commit="tested")))
    with ledger.transaction() as state:
        state["production_gate"]=dict(path=str(gate),sha256=digest_file(gate),producer_commit="tested")
    ledger.reserve("approved",phase="production")
    gate.write_text(json.dumps(dict(passed=False,producer_commit="tested")))
    with pytest.raises(CampaignBudgetError,match="gate artifact"):
        ledger.reserve("corrupt",phase="production")
    gate.unlink()
    with pytest.raises(CampaignBudgetError,match="gate artifact"):
        ledger.reserve("missing",phase="finalisation")
