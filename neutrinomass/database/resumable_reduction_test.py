import gzip
import pytest

from neutrinomass.completions.completions import operator_completions
from neutrinomass.completions.operators import EFF_OPERATORS
from neutrinomass.database.serialization import dumps_completion
from neutrinomass.database.deduplication import deduplicate_completion_jsonl
from neutrinomass.database.checkpoints import WorkBudget, BudgetExpired
from neutrinomass.database.resumable_reduction import reduce_resumable


def test_resumed_reduction_matches_reference_bytes_and_digests(tmp_path):
    source = tmp_path/"raw.jsonl"
    source.write_text("".join(dumps_completion(c)+"\n"
        for c in operator_completions(EFF_OPERATORS["1"])))
    reference = tmp_path/"reference.jsonl"
    expected = deduplicate_completion_jsonl(source,reference)
    directory = tmp_path/"reduction"
    for stop in (1,2,3):
        try:
            reduce_resumable([source],directory,{"seed":0},budget=WorkBudget(fault_after=stop))
        except BudgetExpired:
            pass
    result = reduce_resumable([source],directory,{"seed":0})
    with gzip.open(directory/result["artifacts"]["output"]["path"],"rb") as stream:
        assert stream.read() == reference.read_bytes()
    assert result["progress"]["exact_classes"] == expected["exact_classes"]
    assert result["progress"]["input_completion_digest"] == expected["input_completion_digest"]
    assert result["progress"]["exact_completion_digest"] == expected["exact_completion_digest"]
    assert reduce_resumable([source],directory,{"seed":0}) == result


def test_reducer_rejects_changed_input(tmp_path):
    source = tmp_path/"raw.jsonl"
    source.write_text("")
    reduce_resumable([source],tmp_path/"work",{})
    source.write_text("changed\n")
    with pytest.raises(ValueError,match="Incompatible"):
        reduce_resumable([source],tmp_path/"work",{})
