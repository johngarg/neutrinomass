import json
from types import SimpleNamespace
import pytest

import erratum_submit as controller
from neutrinomass.database.campaign import CampaignLedger
from neutrinomass.database.work_queue import WorkQueue


def configuration(tmp_path):
    config=dict(root=str(tmp_path/"run"),queue=str(tmp_path/"queue.sqlite3"),
        ledger=str(tmp_path/"ledger.json"),phase="validation",campaign_id="test",
        source="/tested/source",interpreter="/pinned/python",seeds=[0,1])
    path=tmp_path/"config.json"
    path.write_text(json.dumps(config))
    return path,config


@pytest.mark.parametrize("limit",[1,3,4])
def test_submission_preserves_reservations_and_respects_ready_dependencies(tmp_path,monkeypatch,limit):
    path,config=configuration(tmp_path)
    config["worker_limit"]=limit
    path.write_text(json.dumps(config))
    queue=WorkQueue(config["queue"])
    for i in range(8):
        queue.add(str(i),{},phase="validation",seed=0)
    queue.add("final",{},phase="validation",seed=0,dependencies=[str(i) for i in range(8)])
    root=tmp_path/"run";root.mkdir();(root/"prepared.json").write_text("{}")
    calls=[]
    def sbatch(command,**kwargs):
        calls.append(command)
        return SimpleNamespace(returncode=0,stdout=str(100+len(calls))+"\n",stderr="")
    monkeypatch.setattr(controller.subprocess,"run",sbatch)
    monkeypatch.setattr(controller,"accounting",lambda ids:{})
    report=controller.tick(path)
    assert len(report["submitted"])==limit
    assert report["accounting"]["charged_core_hours"]==2*limit
    assert not controller.tick(path)["submitted"]
    assert len(calls)==limit
    assert all("--export=ALL,PYTHONHASHSEED=0" in call for call in calls)
    assert "exec /pinned/python -u erratum_worker.py work" in (root/"launch-work.sh").read_text()


def test_saved_continuation_is_charged_but_scientific_failure_blocks(tmp_path,monkeypatch):
    path,config=configuration(tmp_path)
    queue=WorkQueue(config["queue"]);ledger=CampaignLedger(config["ledger"])
    ledger.reserve("first");ledger.submitted("first","100")
    with ledger.transaction() as state:
        state["jobs"]["first"]["queue_config"]=str(path)
    monkeypatch.setattr(controller,"accounting",lambda ids:{"100":dict(status="FAILED",cpus=1,elapsed=300,exit_code="75:0")})
    controller.reconcile(ledger,queue,path)
    assert ledger.summary()["jobs"]["first"]["status"]=="CONTINUED"
    assert ledger.summary()["charged_core_hours"]==300/3600
    ledger.reserve("second");ledger.submitted("second","101")
    with ledger.transaction() as state:
        state["jobs"]["second"]["queue_config"]=str(path)
    monkeypatch.setattr(controller,"accounting",lambda ids:{"101":dict(status="FAILED",cpus=1,elapsed=10,exit_code="1:0")})
    controller.reconcile(ledger,queue,path)
    assert ledger.summary()["submission_block"]["work_id"]=="second"
    with pytest.raises(RuntimeError,match="SUBMISSION_BLOCKED"):
        controller.submit(ledger,config,path,command="work",seed=0)


def test_unbound_submission_never_causes_duplicate_job(tmp_path):
    path,config=configuration(tmp_path)
    ledger=CampaignLedger(config["ledger"]);queue=WorkQueue(config["queue"])
    ledger.reserve("uncertain")
    with ledger.transaction() as state:
        state["jobs"]["uncertain"]["queue_config"]=str(path)
    with pytest.raises(RuntimeError,match="AMBIGUOUS_SUBMISSION"):
        controller.reconcile(ledger,queue,path)


def test_pilot_worker_is_not_a_retry_of_an_unrelated_validation_job(tmp_path,monkeypatch):
    path,config=configuration(tmp_path)
    ledger=CampaignLedger(config["ledger"])
    ledger.reserve("fault-replica")
    ledger.finished("fault-replica",elapsed_seconds=1,allocated_cpus=1,status="NODE_FAIL",failure_class="infrastructure")
    with ledger.transaction() as state:
        state["jobs"]["fault-replica"].update(queue_config=str(path),seed=0,command="recovery-generation")
    monkeypatch.setattr(controller.subprocess,"run",lambda *a,**k:SimpleNamespace(returncode=0,stdout="102\n",stderr=""))
    result=controller.submit(ledger,config,path,command="work",seed=0)
    assert ledger.summary()["jobs"][result["work_id"]]["retry_of"] is None
