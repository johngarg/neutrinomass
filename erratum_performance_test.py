import json
from pathlib import Path

import pytest

from neutrinomass.database.checkpoints import CheckpointStore, WorkBudget, BudgetExpired, digest_file
from neutrinomass.database.resumable_reduction import reduce_resumable
from neutrinomass.database.serialization import dumps_completion
from neutrinomass.completions.operators import EFF_OPERATORS
from neutrinomass.completions.completions import operator_completions
from scripts.erratum_performance_worker import adopt_partial, install, verify_backend
from neutrinomass.database.checkpoints import atomic_json


def backend_fixture(tmp_path):
    from neutrinomass.completions import equivalence
    from neutrinomass.database import global_equivalence, resumable_reduction
    baseline={name:digest_file(module.__file__) for name,module in (
        ('equivalence',equivalence),('global_equivalence',global_equivalence),
        ('resumable_reduction',resumable_reduction))}
    backend=dict(path=equivalence.__file__,sha256=digest_file(equivalence.__file__))
    certificate=tmp_path/'certificate.json'
    atomic_json(certificate,dict(passed=True,backend_sha256=backend['sha256'],baseline_files=baseline,
        evidence=[],full_local_suites=[0,1],full_spartan_suites=[0,1],full_ordered_71p_replays=[0,1]))
    return dict(baseline_files=baseline,module=backend,
                validation=dict(path=str(certificate),sha256=digest_file(certificate)))


def test_installed_backend_records_mixed_source_without_changing_physical_stage(tmp_path,monkeypatch):
    from neutrinomass.database import global_equivalence,pipeline,resumable_reduction
    for module,name in ((global_equivalence,'equivalent_lagrangians'),(pipeline,'reduce_resumable'),
                        (resumable_reduction,'reduce_resumable'),(pipeline,'atomic_json')):
        monkeypatch.setattr(module,name,getattr(module,name))
    spec=backend_fixture(tmp_path)
    install(dict(performance_backend=spec,performance_archive=str(tmp_path/'archive')))
    source=tmp_path/'input.jsonl'
    record=next(operator_completions(EFF_OPERATORS['1']))
    source.write_text((dumps_completion(record)+'\n')*2)
    physical=pipeline.reduce_resumable([source],tmp_path/'physical',dict(source_commit='base'))
    global_=pipeline.reduce_resumable([source],tmp_path/'global',dict(source_commit='base'),relation='global-lagrangian')
    assert 'performance_backend' not in physical['contract']
    assert global_['contract']['source_commit']=='base' and global_['contract']['performance_backend']==spec
    assert physical['artifacts']['output']==global_['artifacts']['output']
    report=tmp_path/'report.json';pipeline.atomic_json(report,dict(stage_manifest=dict(contract='original')))
    assert json.loads(report.read_text())['stage_manifest']['performance_backends']['global_lagrangians']==spec


@pytest.mark.parametrize('key',['full_local_suites','full_spartan_suites','full_ordered_71p_replays'])
def test_backend_requires_each_validation_seed(tmp_path,key):
    spec=backend_fixture(tmp_path);p=Path(spec['validation']['path']);d=json.loads(p.read_text())
    d[key]=[0];atomic_json(p,d);spec['validation']['sha256']=digest_file(p)
    with pytest.raises(ValueError,match='PERFORMANCE_VALIDATION_INCOMPLETE'):
        verify_backend(spec)


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
