import pytest
from neutrinomass.database.work_queue import WorkQueue


def test_dependency_queue_does_not_duplicate_publication(tmp_path):
    queue=WorkQueue(tmp_path/"queue.sqlite3")
    queue.add("generate",{},phase="pilot",seed=0)
    queue.add("finish",{},phase="pilot",seed=0,dependencies=["generate"])
    assert queue.claim("job1",phase="pilot",seed=0)["id"] == "generate"
    assert queue.claim("job2",phase="pilot",seed=0) is None
    queue.finish("generate","job1",{"sha256":"x"})
    queue.finish("generate","job1",{"sha256":"x"})
    assert queue.claim("job2",phase="pilot",seed=0)["id"] == "finish"
    with pytest.raises(ValueError):
        queue.finish("finish","wrong-job",{})


def test_repeated_no_progress_stops_the_unit(tmp_path):
    queue=WorkQueue(tmp_path/"queue.sqlite3")
    queue.add("unit",{},phase="pilot",seed=0)
    for expected in ("pending","pending","blocked"):
        queue.claim("job",phase="pilot",seed=0)
        assert queue.continuation("unit","job",{"cursor":4}) == expected
    assert queue.claim("another",phase="pilot",seed=0) is None
    assert queue.summary()["blocked"][0]["diagnostic"] == "REPEATED_NO_PROGRESS"


def test_infrastructure_retries_are_bounded_and_scientific_failures_block(tmp_path):
    queue=WorkQueue(tmp_path/"queue.sqlite3")
    queue.add("unit",{},phase="pilot",seed=0)
    for i in range(3):
        queue.claim(str(i),phase="pilot",seed=0)
        queue.recover_infrastructure(str(i),"NODE_FAIL")
    assert queue.claim("again",phase="pilot",seed=0) is None
    with pytest.raises(ValueError):
        queue.recover_infrastructure("0","FAILED")
