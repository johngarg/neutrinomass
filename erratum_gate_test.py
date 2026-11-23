from copy import deepcopy
import pytest
from scripts.validate_71p_gate import report_gate,recovery_gate,GateFailure
from neutrinomass.database.checkpoints import CheckpointStore,atomic_json,digest_file


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
