from types import SimpleNamespace

import pytest

from census_derivative_operator import validate_physics, verify_round_trip
from neutrinomass.completions.completions import operator_completions
from neutrinomass.completions.fingerprints import completion_digest
from neutrinomass.completions.operators import EFF_OPERATORS
from neutrinomass.database import write_completion_jsonl
from neutrinomass.tensormethod import H, eps


def test_validate_physics_rejects_vanishing_uv_interaction():
    vanishing = H("i0") * H("i1") * eps("-i0 -i1")
    completion = SimpleNamespace(terms=[vanishing], derivative_routes=())

    with pytest.raises(ValueError, match="vanishing UV interaction"):
        validate_physics([completion])


def test_verify_round_trip_streams_completion_digest(tmp_path):
    completion = next(operator_completions(EFF_OPERATORS["1"]))
    path = tmp_path / "completion.jsonl"
    write_completion_jsonl(path, [completion])

    assert verify_round_trip(path, [completion]) == completion_digest([completion])
