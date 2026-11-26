import json
from pathlib import Path

import pytest

from neutrinomass.database.checkpoints import CheckpointStore, WorkBudget, BudgetExpired, digest_file
from neutrinomass.database.resumable_reduction import reduce_resumable
from neutrinomass.database.serialization import dumps_completion
from neutrinomass.completions.operators import EFF_OPERATORS
from neutrinomass.completions.completions import operator_completions
from scripts.erratum_performance_worker import adopt_partial


def test_performance_transition_preserves_cursor_and_ordered_result(tmp_path):
    records=list(operator_completions(EFF_OPERATORS['1']))
    source=tmp_path/'input.jsonl'
    source.write_text('\n'.join(dumps_completion(c) for c in records*2)+'\n')
    contract=dict(source_commit='original')
    full=reduce_resumable([source],tmp_path/'full',contract,relation='global-lagrangian')
    with pytest.raises(BudgetExpired):
        reduce_resumable([source],tmp_path/'partial',contract,relation='global-lagrangian',
                         budget=WorkBudget(fault_after=4))
    path=tmp_path/'partial/checkpoint.json';old=json.loads(path.read_text())
    old_hash=digest_file(path);spec=dict(validation=dict(path='certificate',sha256='tested'))
    new_identity=dict(old['contract'],performance_backend=spec)
    adopt_partial(path.parent,old['contract'],new_identity,tmp_path/'archive',spec)
    changed=json.loads(path.read_text())
    assert changed['progress']==old['progress'] and changed['artifacts']==old['artifacts']
    assert digest_file(tmp_path/'archive'/old_hash/'checkpoint.json')==old_hash
    # Repeat after an interrupted launch: migration is idempotent.
    adopt_partial(path.parent,old['contract'],new_identity,tmp_path/'archive',spec)
    continued=reduce_resumable([source],path.parent,dict(contract,performance_backend=spec),
                               relation='global-lagrangian')
    assert continued['progress']==full['progress']
    assert continued['artifacts']['output']==full['artifacts']['output']
    with CheckpointStore(tmp_path/'archive'/old_hash,old['contract']) as archived:
        assert archived.state['progress']==old['progress']


@pytest.mark.parametrize('damage',['artifact','identity','complete'])
def test_performance_transition_rejects_incompatible_or_damaged_stage(tmp_path,damage):
    identity=dict(source_commit='original',relation='global-lagrangian')
    with CheckpointStore(tmp_path/'stage',identity) as store:
        artifact=store.publish_lines(['record'])
        store.save(dict(processed=1),artifacts=dict(output=artifact),complete=damage=='complete')
    spec=dict(validation={});new=dict(identity,performance_backend=spec)
    if damage=='artifact':
        (tmp_path/'stage'/artifact['path']).write_text('corrupt')
    elif damage=='identity':
        new['source_commit']='different'
    with pytest.raises(ValueError):
        adopt_partial(tmp_path/'stage',identity,new,tmp_path/'archive',spec)
