from copy import deepcopy
import pytest
from scripts.validate_71p_gate import report_gate,recovery_gate,GateFailure
from neutrinomass.database.checkpoints import CheckpointStore,atomic_json,digest_file


def test_identical_seed_bytes_prove_agreement_without_reducing_again(tmp_path,monkeypatch):
    from scripts import validate_71p_gate as gate
    reports=[]
    for seed in (0,1):
        p=tmp_path/f"seed-{seed}.jsonl";p.write_text('same complete serialized record\n')
        ref=dict(path=str(p),sha256=digest_file(p))
        reports.append(dict(exact_classes=dict(all=1),global_lagrangians=1,
            artifacts=dict(exact_classes=ref,global_lagrangians=ref)))
    def unexpected(*args,**kwargs):
        raise AssertionError("Identical bytes need no repeated scientific computation")
    monkeypatch.setattr(gate,"historical_coverage",unexpected)
    monkeypatch.setattr(gate,"reduce_resumable",unexpected)
    exact,union=gate.compare_seed_artifacts(reports,tmp_path,{},None)
    assert exact[0]["counts"]==dict(reproduced=1)
    assert union["method"]=="identical-serialized-bytes" and union["exact_classes"]==1
    # A stale checksum must never certify altered/truncated bytes.
    (tmp_path/'seed-1.jsonl').write_text('changed')
    with pytest.raises(GateFailure,match='GATE_ARTIFACT_HASH_MISMATCH'):
        gate.compare_seed_artifacts(reports,tmp_path,{},None)


def test_different_seed_bytes_keep_both_exact_coverage_checks_and_union(tmp_path,monkeypatch):
    from scripts import validate_71p_gate as gate
    reports=[];calls=[]
    for seed in (0,1):
        p=tmp_path/f"seed-{seed}.jsonl";p.write_text(str(seed))
        ref=dict(path=str(p),sha256=digest_file(p))
        reports.append(dict(exact_classes=dict(all=1),global_lagrangians=1,
            artifacts=dict(exact_classes=ref,global_lagrangians=ref)))
    def coverage(*args,**kwargs):
        calls.append('coverage');return dict(progress=dict(counts=dict(reproduced=1)))
    def union(*args,**kwargs):
        calls.append('union');return dict(progress=dict(exact_classes=1))
    monkeypatch.setattr(gate,'historical_coverage',coverage)
    monkeypatch.setattr(gate,'reduce_resumable',union)
    monkeypatch.setattr(gate,'ordered_seed_agreement',lambda *args:dict(agreement=False))
    gate.compare_seed_artifacts(reports,tmp_path,{},None)
    assert calls==['coverage','coverage','union']


@pytest.mark.parametrize('relation',['physical','global-lagrangian'])
def test_ordered_seed_proof_resumes_and_rejects_reordered_pair(tmp_path,relation):
    from scripts.validate_71p_gate import ordered_seed_agreement
    from neutrinomass.database.checkpoints import WorkBudget,BudgetExpired
    from neutrinomass.completions.completions import operator_completions
    from neutrinomass.completions.operators import EFF_OPERATORS
    from neutrinomass.database.serialization import dumps_completion
    records=list(operator_completions(EFF_OPERATORS['1']))
    refs=[]
    for seed in (0,1):
        p=tmp_path/f'seed-{seed}.jsonl';p.write_text('\n'.join(map(dumps_completion,records))+'\n')
        refs.append(dict(path=str(p),sha256=digest_file(p)))
    with pytest.raises(BudgetExpired):
        ordered_seed_agreement(refs,tmp_path/'proof',{},relation,WorkBudget(fault_after=1))
    result=ordered_seed_agreement(refs,tmp_path/'proof',{},relation,WorkBudget())
    assert result==dict(agreement=True,compared_pairs=len(records),comparisons=2*len(records))
    p=tmp_path/'seed-1.jsonl';p.write_text('\n'.join(map(dumps_completion,reversed(records)))+'\n')
    changed=[refs[0],dict(path=str(p),sha256=digest_file(p))]
    with pytest.raises(ValueError,match='Incompatible checkpoint'):
        ordered_seed_agreement(changed,tmp_path/'proof',{},relation,WorkBudget())
    result=ordered_seed_agreement(changed,tmp_path/'different-order',{},relation,WorkBudget())
    assert not result['agreement']


def report():
    return dict(operator="71p",hash_seed="0",source_commit="tested",records=dict(generator=758208),
        stage_manifest=dict(generation_units=3168,topologies=[str(i) for i in range(264)]),
        historical=dict(missing=0,amplitude_unsupported=[],valid_classes=2,reproduced=2),
        exact_classes=dict(all=3))


@pytest.mark.parametrize("mutation",[
    lambda r:r["stage_manifest"]["topologies"].pop(),
    lambda r:r["records"].update(generator=758207),
    lambda r:r["historical"].update(missing=1),
    lambda r:r["historical"].update(amplitude_unsupported=["unresolved"]),
    lambda r:r.update(hash_seed="1"),
    lambda r:r.update(source_commit="untested"),
])
def test_incomplete_or_unexplained_pilot_cannot_open_production(mutation):
    candidate=report()
    report_gate(candidate,seed=0,producer="tested",unit_count=3168)
    mutation(candidate)
    with pytest.raises(GateFailure):
        report_gate(candidate,seed=0,producer="tested",unit_count=3168)


@pytest.mark.parametrize("damage",["wrong-operator","corrupt-chunk"])
def test_recovery_gate_checks_the_actual_archived_fault_snapshot(tmp_path,damage):
    archive=tmp_path/"interruptions/job"
    with CheckpointStore(archive,dict(source_commit="tested",operator="1" if damage=="wrong-operator" else "71p")) as store:
        chunk=store.publish_lines(["durable"])
        store.save(dict(next_assignment=5),artifacts={"chunk":chunk})
        progress=store.state["progress"]
    checkpoint=archive/"checkpoint.json"
    atomic_json(tmp_path/"interruptions/job.json",dict(kind="generation",unit=0,injected=True,
        method="SIGKILL",child_returncode=-9,checkpoints=[dict(path=str(checkpoint),
        sha256=digest_file(checkpoint),progress=progress)]))
    atomic_json(tmp_path/"recovery-evidence.json",dict(complete=True,producer=dict(source_commit="tested"),
        job_id="job",kind="generation",unit=0,ordered_exact_agreement=True))
    if damage=="corrupt-chunk":
        (archive/chunk["path"]).write_bytes(b"truncated")
    with pytest.raises(ValueError,match="GATE_FAULT_OPERATOR_MISMATCH|Damaged checkpoint artifact"):
        recovery_gate([tmp_path],"tested")
