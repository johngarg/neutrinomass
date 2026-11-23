import pytest

from neutrinomass.database.filtering import (
    filter_published_registry, filter_species_registry, multiset_contains,
    published_model_key,
)
from neutrinomass.utils.functions import conjugate_field


def test_published_identity_preserves_conjugation_and_repeated_factors():
    field = "S,10,0,1/3,1"
    assert published_model_key([field]) == published_model_key([conjugate_field(field)])
    assert published_model_key([field,field]) != published_model_key([field])


def test_published_mass_witness_respects_equal_scale_groups():
    a, b = "S,00,0,7,0", "S,00,0,8,0"
    result = filter_published_registry(
        {"high": {(a,): ()}, "equal": {(a,b): ()}, "low": {(a,b): ()}},
        {"high": 3.0, "equal": 3.0, "low": 1.0})
    assert result["survivors"]["equal"]
    assert not result["survivors"]["low"]
    assert result["witnesses"][0]["code"] == "higher_scale_published_model"


def test_two_copies_cannot_exclude_one_copy_and_incomplete_search_raises():
    assert multiset_contains(("S","S","F"),("S","S"))
    assert not multiset_contains(("S","F"),("S","S"))
    with pytest.raises(ValueError, match="absent"):
        filter_species_registry({"x": {("S","F"): ()}}, {"x": 1},
            mechanism_search=lambda *args: {"required_species": ["S","S"], "verified": True},
            loop_search=lambda *args: None)
    def incomplete(*args):
        raise RuntimeError("Incomplete candidate mechanism search")
    with pytest.raises(RuntimeError, match="Incomplete"):
        filter_species_registry({"x": {("S",): ()}}, {"x": 1},
            mechanism_search=incomplete, loop_search=lambda *args: None)
