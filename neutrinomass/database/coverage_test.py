import json
import pytest

from neutrinomass.database.deduplication_test import relabelled_scalar_completion
from neutrinomass.database.serialization import dumps_completion,loads_completion
from neutrinomass.database.checkpoints import WorkBudget,BudgetExpired
from neutrinomass.database.coverage import historical_coverage
from neutrinomass.database.resumable_reduction import lines_from_paths,reduce_resumable
from neutrinomass.database.pipeline import artifact_path
from neutrinomass.tensormethod import H,eps


def bucket_pair():
    first=relabelled_scalar_completion("a","i0","i1")
    second=loads_completion(dumps_completion(first))
    second.terms=[first.terms[0]*H("i3")*H.conj("i4")*eps("-i3 -i4")]
    return first,second


def test_coverage_resumes_inside_a_bucket_without_repeating_comparisons(tmp_path):
    first,second=bucket_pair()
    exact=tmp_path/"exact.jsonl"
    historical=tmp_path/"historical.jsonl"
    exact.write_text(dumps_completion(first)+"\n"+dumps_completion(second)+"\n")
    historical.write_text(dumps_completion(second)+"\n")
    full=historical_coverage([exact],[historical],tmp_path/"full",{})
    with pytest.raises(BudgetExpired):
        historical_coverage([exact],[historical],tmp_path/"resumed",{},budget=WorkBudget(fault_after=3))
    saved=json.loads((tmp_path/"resumed/checkpoint.json").read_text())
    assert saved["progress"]["pending"]["cursor"] == 1
    resumed=historical_coverage([exact],[historical],tmp_path/"resumed",{})
    assert resumed["progress"]["comparisons"] == full["progress"]["comparisons"] == 2
    assert resumed["progress"]["counts"] == {"reproduced":1}
    assert list(lines_from_paths([artifact_path(tmp_path/"full",full)])) == list(lines_from_paths([artifact_path(tmp_path/"resumed",resumed)]))


def test_reduction_can_stop_after_a_failed_exact_comparison(tmp_path,monkeypatch):
    import neutrinomass.database.resumable_reduction as reducer
    first,second=bucket_pair()
    # A colliding search index must not change the equality decision.
    monkeypatch.setattr(reducer,"completion_fingerprint",lambda c:(0,0,0,0,"collision"))
    source=tmp_path/"records.jsonl"
    source.write_text("\n".join(map(dumps_completion,[first,second,second]))+"\n")
    full=reduce_resumable([source],tmp_path/"full",{})
    with pytest.raises(BudgetExpired):
        reduce_resumable([source],tmp_path/"resumed",{},budget=WorkBudget(fault_after=2))
    saved=json.loads((tmp_path/"resumed/checkpoint.json").read_text())
    assert saved["progress"]["pending"]["ordinal"] == 2
    assert saved["progress"]["pending"]["cursor"] == [0,1]
    resumed=reduce_resumable([source],tmp_path/"resumed",{})
    assert resumed["progress"]["exact_classes"] == 2
    assert resumed["progress"]["comparisons"] == full["progress"]["comparisons"]
    assert list(lines_from_paths([artifact_path(tmp_path/"full",full)])) == list(lines_from_paths([artifact_path(tmp_path/"resumed",resumed)]))
