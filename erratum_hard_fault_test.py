import json
import signal
import sys
import pytest
from scripts import erratum_hard_fault


@pytest.mark.parametrize("kind",["generation","reduction"])
def test_replica_is_killed_and_only_persistent_artifacts_reach_new_process(tmp_path,monkeypatch,kind):
    driver=tmp_path/"erratum_hard_fault.py"
    driver.write_text("test driver identity\n")
    child=tmp_path/"erratum_recovery.py"
    child.write_text('''import json,sys,time,os
from pathlib import Path
replica=Path(sys.argv[4]); kind=sys.argv[1]
stage=replica/'units/00000000/generation' if kind=='generation' else replica
if '--checkpoint-seconds' in sys.argv:
    stage.mkdir(parents=True)
    (stage/'piece').write_text('durable')
    progress=dict(next_assignment=5) if kind=='generation' else dict(processed=2,pending=dict(ordinal=3,cursor=[0,1]))
    (stage/'checkpoint.json').write_text(json.dumps(dict(complete=False,progress=progress,artifacts=dict(piece=dict(path='piece')))))
    while True: time.sleep(1)
else:
    assert (stage/'piece').read_text()=='durable'
    (replica/'recovered-process.json').write_text(json.dumps(dict(pid=os.getpid())))
''')
    replica=tmp_path/"replica"
    monkeypatch.setattr(erratum_hard_fault,"__file__",str(driver))
    monkeypatch.setenv("SLURM_JOB_ID","test-job")
    monkeypatch.setattr(sys,"argv",[str(driver),kind,str(tmp_path/"config"),str(tmp_path/"reference"),str(replica),"--unit","0"])
    previous={sig:signal.getsignal(sig) for sig in (signal.SIGUSR1,signal.SIGTERM)}
    try:
        assert erratum_hard_fault.main()==0
    finally:
        for sig,handler in previous.items():signal.signal(sig,handler)
    evidence=json.loads((replica/"interruptions/test-job.json").read_text())
    assert evidence["method"]=="SIGKILL" and evidence["child_returncode"]==-9
    assert evidence["in_bucket_reduction"]==(kind=="reduction")
    assert (replica/"interruptions/test-job/piece").read_text()=="durable"
    assert (replica/"recovered-process.json").exists()
    # A later Slurm allocation resumes the same replica without injecting a
    # second fault or overwriting the original durable fault certificate.
    old=(replica/"interruptions/test-job.json").read_bytes()
    monkeypatch.setenv("SLURM_JOB_ID","next-job")
    try:
        assert erratum_hard_fault.main()==0
    finally:
        for sig,handler in previous.items():signal.signal(sig,handler)
    assert (replica/"interruptions/test-job.json").read_bytes()==old
    assert not (replica/"interruptions/next-job.json").exists()
