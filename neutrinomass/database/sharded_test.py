import json
import pytest

from neutrinomass.database.checkpoints import WorkBudget,BudgetExpired
from neutrinomass.database.coverage import historical_coverage
from neutrinomass.database.sharded import sharded_map,sharded_coverage
from neutrinomass.database.stages import map_records,amplitude_record
from neutrinomass.database.pipeline import artifact_path
from neutrinomass.database.resumable_reduction import lines_from_paths
from neutrinomass.database.serialization import dumps_completion
from neutrinomass.completions.completions import operator_completions
from neutrinomass.completions.operators import EFF_OPERATORS


def source_records(tmp_path):
    path=tmp_path/"records.jsonl"
    records=list(operator_completions(EFF_OPERATORS["1"]))
    path.write_text("\n".join(map(dumps_completion,records))+"\n")
    return path,len(records)


def test_out_of_order_amplitude_shards_merge_in_original_order(tmp_path):
    path,count=source_records(tmp_path)
    full=map_records([path],tmp_path/"full",{},amplitude_record,algorithm="amplitudes-v1")
    sharded_map([path],count,tmp_path/"sharded",{},amplitude_record,algorithm="amplitudes-v1",
                shard_size=2,selected_unit=3)
    with pytest.raises(BudgetExpired):
        sharded_map([path],count,tmp_path/"sharded",{},amplitude_record,algorithm="amplitudes-v1",
                    shard_size=2,budget=WorkBudget(fault_after=1))
    recovered=sharded_map([path],count,tmp_path/"sharded",{},amplitude_record,algorithm="amplitudes-v1",shard_size=2)
    assert full["progress"]==recovered["progress"]
    for key in ("output","evidence"):
        assert list(lines_from_paths([artifact_path(tmp_path/"full",full,key)])) == list(
            lines_from_paths([artifact_path(tmp_path/"sharded",recovered,key)]))


def test_coverage_shards_share_immutable_index_and_recover_exact_matches(tmp_path):
    path,count=source_records(tmp_path)
    full=historical_coverage([path],[path],tmp_path/"full",{})
    sharded_coverage([path],[path],count,tmp_path/"sharded",{},shard_size=2,selected_unit=3)
    with pytest.raises(BudgetExpired):
        sharded_coverage([path],[path],count,tmp_path/"sharded",{},shard_size=2,
                         budget=WorkBudget(fault_after=1))
    result=sharded_coverage([path],[path],count,tmp_path/"sharded",{},shard_size=2)
    assert result["progress"]["counts"]==full["progress"]["counts"]=={"reproduced":count}
    assert result["progress"]["comparisons"]==full["progress"]["comparisons"]
    assert list(lines_from_paths([artifact_path(tmp_path/"full",full)])) == list(
        lines_from_paths([artifact_path(tmp_path/"sharded",result)]))
    index=json.loads((tmp_path/"sharded/index/checkpoint.json").read_text())
    for checkpoint in (tmp_path/"sharded/shards").glob("*/checkpoint.json"):
        assert json.loads(checkpoint.read_text())["contract"]["shared_index"][0]["sha256"]==index["artifacts"]["database"]["sha256"]
