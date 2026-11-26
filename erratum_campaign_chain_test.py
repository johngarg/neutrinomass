import json
from types import SimpleNamespace
import pytest

from scripts import run_erratum_campaign as chain_driver
from neutrinomass.database.campaign import CampaignLedger,CampaignBudgetError
from neutrinomass.database.checkpoints import atomic_json


def setup(tmp_path):
    pilot=tmp_path/"pilot.json"
    atomic_json(pilot,dict(interpreter="/pinned/python",source="/tested/source"))
    ledger=CampaignLedger(tmp_path/"ledger.json")
    chain=dict(id="test-chain",root=str(tmp_path/"chain"))
    spec=dict(phase="validation",source="/tested/source",argv=["/tested/script.py","literal $value"],
              progress=[str(tmp_path/"progress")],evidence=[str(tmp_path/"proof.json")])
    return ledger,pilot,chain,spec


def test_chain_submission_is_bounded_and_restart_does_not_duplicate(tmp_path,monkeypatch):
    ledger,pilot,chain,spec=setup(tmp_path);calls=[]
    def submit(command,**kwargs):
        calls.append(command);return SimpleNamespace(returncode=0,stdout="10\n",stderr="")
    monkeypatch.setattr(chain_driver.subprocess,"run",submit)
    assert chain_driver.external(ledger,pilot,chain,"fault",spec)=="running"
    assert chain_driver.external(ledger,pilot,chain,"fault",spec)=="running"
    assert len(calls)==1
    assert ledger.summary()["charged_core_hours"]==2
    batch=(tmp_path/"chain/fault-0000.sh").read_text()
    assert "'literal $value'" in batch
    ledger.finished("test-chain/fault/0000",elapsed_seconds=10,allocated_cpus=1,status="COMPLETED")
    with pytest.raises(ValueError,match="WITHOUT_EVIDENCE"):
        chain_driver.external(ledger,pilot,chain,"fault",spec)
    atomic_json(tmp_path/"proof.json",dict(complete=True))
    assert chain_driver.external(ledger,pilot,chain,"fault",spec)=="complete"
    with pytest.raises(ValueError,match="COMMAND_CHANGED"):
        chain_driver.external(ledger,pilot,chain,"fault",dict(spec,argv=["different"]))


def test_chain_stops_repeated_continuations_without_durable_progress(tmp_path,monkeypatch):
    ledger,pilot,chain,spec=setup(tmp_path)
    monkeypatch.setattr(chain_driver.subprocess,"run",lambda *a,**k:SimpleNamespace(returncode=0,stdout="10\n",stderr=""))
    for attempt in range(3):
        assert chain_driver.external(ledger,pilot,chain,"reduce",spec)=="running"
        ledger.finished(f"test-chain/reduce/{attempt:04d}",elapsed_seconds=1,allocated_cpus=1,status="CONTINUED")
    with pytest.raises(ValueError,match="REPEATED_NO_PROGRESS"):
        chain_driver.external(ledger,pilot,chain,"reduce",spec)
    assert len(ledger.summary()["jobs"])==3
    assert ledger.summary()["submission_block"]["code"]=="CHAIN_REPEATED_NO_PROGRESS"


def test_chain_never_submits_production_before_gate_or_over_concurrency(tmp_path,monkeypatch):
    ledger,pilot,chain,spec=setup(tmp_path)
    with pytest.raises(CampaignBudgetError,match="gate"):
        chain_driver.external(ledger,pilot,chain,"production",dict(spec,phase="production"))
    for ordinal in range(4):ledger.reserve(str(ordinal))
    assert chain_driver.external(ledger,pilot,chain,"fault",spec)=="awaiting-slot"
    assert len(ledger.summary()["jobs"])==4


def test_chain_will_not_resubmit_an_ambiguous_unbound_reservation(tmp_path):
    ledger,pilot,chain,spec=setup(tmp_path)
    ledger.reserve("test-chain/fault/0000")
    with ledger.transaction() as state:
        state["jobs"]["test-chain/fault/0000"].update(chain_stage="fault",chain_id=chain["id"],
            chain_identity=chain_driver.identity_hash(spec))
    # The controller's reconcile raises on this unbound job; external itself
    # also retains the reservation and cannot send a duplicate sbatch.
    assert chain_driver.external(ledger,pilot,chain,"fault",spec)=="running"
    assert len(ledger.summary()["jobs"])==1


@pytest.mark.parametrize('gate_override',[False,True])
def test_active_production_allocations_do_not_block_their_own_accounting_tick(tmp_path,monkeypatch,gate_override):
    ledger,pilot,chain,spec=setup(tmp_path)
    config=dict(root=str(tmp_path/"pilot"),ledger=str(ledger.path),queue=str(tmp_path/"pilot-queue.sqlite3"),
                source="/producer",source_commit="producer",worker_limit=4)
    atomic_json(pilot,config)
    def ref(path):return dict(path=str(path),sha256=chain_driver.digest_file(path))
    certificates=[]
    for seed in (0,1):
        path=tmp_path/f"test-{seed}.json"
        atomic_json(path,dict(full_suite_passed=True,environment=dict(source_commit="finaliser",dirty_worktree=False,hash_seed=str(seed))))
        certificates.append(ref(path))
    generation=tmp_path/"generation"
    atomic_json(generation/"recovery-evidence.json",dict(complete=True,ordered_exact_agreement=True))
    atomic_json(tmp_path/"pilot/seed-0/71p/structural/checkpoint.json",dict(complete=True))
    atomic_json(tmp_path/"pilot/seed-0/71p/report.json",{})
    baseline=tmp_path/"baseline.json";atomic_json(baseline,{})
    chain.update(pilot_config=str(pilot),scripts=str(tmp_path/"scripts"),orchestrator_sources=[],
        finaliser_commit="finaliser",finaliser_validation=certificates,generation_recovery=str(generation),
        baseline=ref(baseline),baseline_models=ref(baseline),production_root=str(tmp_path/"production"))
    ledger.reserve("production-worker")
    with ledger.transaction() as state:
        state["jobs"]["production-worker"].update(command="work",queue_config=str(tmp_path/"production/production-config.json"))
    monkeypatch.setattr(chain_driver,"reconcile",lambda *a:None)
    submissions=[]
    monkeypatch.setattr(chain_driver,"external",lambda *a:submissions.append(a) or "complete")
    if gate_override:
        script=tmp_path/'certified-gate.py';script.write_text('# verified gate\n')
        chain['gate_script']=ref(script)
    monkeypatch.setattr(chain_driver,"tick",lambda *a:dict(status="complete"))
    called=[]
    monkeypatch.setattr(chain_driver,"queue_tick",lambda path:called.append(path) or dict(status="running"))
    result=chain_driver.step(chain)
    assert result["status"]=="running" and result["phase"]=="production"
    assert called==[tmp_path/"production/production-config.json"]
    gate_spec=next(args[-1] for args in submissions if args[-2]=='gate')
    assert gate_spec['argv'][0]==str(tmp_path/'certified-gate.py' if gate_override else tmp_path/'scripts/validate_71p_gate.py')
    if gate_override:
        script.write_text('# changed\n')
        with pytest.raises(ValueError,match='CHAIN_REFERENCE_CHANGED'):
            chain_driver.step(chain)
