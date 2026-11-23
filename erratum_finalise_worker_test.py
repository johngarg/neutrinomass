import json
from pathlib import Path
import pytest

from erratum_finalise_worker import prepare,expand_species_queue
from neutrinomass.database.checkpoints import CheckpointStore,WorkBudget,BudgetExpired,canonical_json
from neutrinomass.database.work_queue import WorkQueue


def test_filter_queue_expansion_is_recoverable_and_final_publication_waits(tmp_path):
    config=dict(root=str(tmp_path),queue=str(tmp_path/"queue.sqlite3"),_sha256="control")
    prepare(config)
    queue=WorkQueue(config["queue"])
    selected=queue.claim("verify-worker",phase="finalisation",seed=0)
    assert selected["id"]=="verify"
    queue.finish(selected["id"],"verify-worker",{})
    selected=queue.claim("published-worker",phase="finalisation",seed=0)
    assert selected["id"]=="published"
    selected=queue.claim("index-worker",phase="finalisation",seed=0)
    assert selected["id"]=="species-index"
    with CheckpointStore(tmp_path/"species-index",{}) as store:
        output=store.publish_lines(canonical_json(dict(identifier=str(i),species="[]",scale_tev=i)) for i in range(3))
        store.save({},artifacts={"work":output},complete=True)
        state=store.state
    with pytest.raises(BudgetExpired):
        expand_species_queue(config,state,{},budget=WorkBudget(fault_after=2))
    assert not any(row["status"]=="blocked" for row in queue.summary()["counts"])
    complete=expand_species_queue(config,state,{},budget=WorkBudget())
    assert complete["progress"]["added"]==3
    assert expand_species_queue(config,state,{},budget=WorkBudget())==complete
    with queue.connection() as db:
        assert db.execute("SELECT COUNT(*) FROM units WHERE id LIKE 'species/%'").fetchone()[0]==3
        parents={row[0] for row in db.execute("SELECT dependency FROM dependencies WHERE unit='species-collect'")}
        assert parents=={"species-index","species/0","species/1","species/2"}
        report=db.execute("SELECT status FROM units WHERE id='report'").fetchone()[0]
        assert report=="pending"
        ready=db.execute("""SELECT COUNT(*) FROM units u WHERE u.id='report' AND NOT EXISTS
            (SELECT 1 FROM dependencies d JOIN units p ON p.id=d.dependency
             WHERE d.unit=u.id AND p.status!='complete')""").fetchone()[0]
        assert ready==0
