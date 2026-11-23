import json
import erratum_worker as worker
from neutrinomass.database.work_queue import WorkQueue


def test_queue_worker_resumes_when_another_worker_publishes_the_plan(tmp_path,monkeypatch):
    config=dict(queue=str(tmp_path/"queue.sqlite3"),phase="validation",root=str(tmp_path/"run"),
                historical=str(tmp_path),assignment_size=2)
    queue=WorkQueue(config["queue"])
    queue.add("first",dict(operator="1",unit=0,stage="local-reduction"),phase="validation",seed=0)
    provenance={key:None for key in ("source_commit","dirty_worktree","working_tree_digest",
        "lockfile_sha256","dependencies","python_version","hash_seed")}
    monkeypatch.setenv("SLURM_JOB_ID","test-worker")
    monkeypatch.setenv("SLURM_TMPDIR",str(tmp_path))
    monkeypatch.setenv("PYTHONHASHSEED","0")
    def raced(operator,historical,output,contract,**kwargs):
        output.mkdir(parents=True,exist_ok=True)
        (output/"plan.json").write_text(json.dumps(dict(created_by="another worker")))
        assert kwargs["resume"] is True
        return {}
    monkeypatch.setattr(worker,"run_operator",raced)
    assert worker.work(config,provenance)==0
    assert queue.summary()["counts"][0]["status"]=="complete"
