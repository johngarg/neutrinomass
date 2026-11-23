import json
import sqlite3
import pytest

from neutrinomass.database.checkpoints import CheckpointStore, WorkBudget, BudgetExpired


def test_sealed_database_survives_unpublished_mutation(tmp_path):
    contract = dict(seed=0, source="fixed", input_sha256="abc")
    with CheckpointStore(tmp_path,contract) as store:
        with store.database() as connection:
            connection.execute("CREATE TABLE results (value INTEGER)")
            connection.execute("INSERT INTO results VALUES (1)")
            store.snapshot(connection,dict(cursor=1))
            connection.execute("INSERT INTO results VALUES (2)")
            connection.commit()
    with CheckpointStore(tmp_path,contract) as store:
        assert store.state["progress"]["cursor"] == 1
        with store.database() as connection:
            assert connection.execute("SELECT * FROM results").fetchall() == [(1,)]


def test_resume_rejects_changed_input_and_damaged_chunk(tmp_path):
    with CheckpointStore(tmp_path,{"input":"a"}) as store:
        artifact = store.publish_lines(['{"a":1}'])
        store.save({"cursor":1},artifacts={"chunk":artifact})
    with pytest.raises(ValueError,match="Incompatible"):
        with CheckpointStore(tmp_path,{"input":"b"}):
            pass
    (tmp_path/artifact["path"]).write_bytes(b"truncated")
    with pytest.raises(ValueError,match="Damaged"):
        with CheckpointStore(tmp_path,{"input":"a"}):
            pass


def test_concurrent_writer_is_rejected_and_orphan_ignored(tmp_path):
    with CheckpointStore(tmp_path,{}) as first:
        first.publish_lines(["orphan"])
        with pytest.raises(BlockingIOError):
            with CheckpointStore(tmp_path,{}):
                pass
        first.save({"cursor":0})
    with CheckpointStore(tmp_path,{}) as second:
        assert not second.state["artifacts"]


def test_time_budget_publishes_before_continuation(tmp_path):
    with CheckpointStore(tmp_path,{}) as store:
        budget = WorkBudget(fault_after=1)
        with pytest.raises(BudgetExpired):
            budget.boundary(lambda: store.save({"cursor":7}))
    with CheckpointStore(tmp_path,{}) as store:
        assert store.state["progress"]["cursor"] == 7


def test_complete_checkpoint_and_compression_are_deterministic(tmp_path):
    with CheckpointStore(tmp_path,{}) as store:
        first = store.publish_lines(["a","b"])
        second = store.publish_lines(["a\n","b\n"])
        assert first == second
        store.save({},artifacts={"output":first},complete=True)
    with CheckpointStore(tmp_path,{}) as store:
        assert store.state["complete"]
        assert store.state["artifacts"]["output"]["records"] == 2
