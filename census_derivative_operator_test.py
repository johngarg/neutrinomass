from types import SimpleNamespace

import pytest

from census_derivative_operator import validate_physics
from neutrinomass.tensormethod import H, eps


def test_validate_physics_rejects_vanishing_uv_interaction():
    vanishing = H("i0") * H("i1") * eps("-i0 -i1")
    completion = SimpleNamespace(terms=[vanishing], derivative_routes=())

    with pytest.raises(ValueError, match="vanishing UV interaction"):
        validate_physics([completion])
