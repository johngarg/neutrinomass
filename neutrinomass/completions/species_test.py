from pathlib import Path
from collections import Counter
import pytest

from neutrinomass.completions.species import candidate_completion_stream
from neutrinomass.completions.operators import EFF_OPERATORS
from neutrinomass.completions.amplitudes import audit_amplitude_symmetrisation
from neutrinomass.database.serialization import iter_completion_jsonl, dumps_completion, loads_completion
from neutrinomass.database.global_equivalence import species_descriptors
from neutrinomass.database.filtering import filter_species_registry

FIXTURES=Path(__file__).resolve().parents[2]/"tests/fixtures/copy-mechanisms"


@pytest.mark.parametrize("case",range(3))
def test_saved_copy_dependent_mechanism_is_generated_before_vertex_rejection(case):
    witness=next(iter_completion_jsonl(FIXTURES/f"{case}.jsonl"))
    assert audit_amplitude_symmetrisation(witness).status == "nonzero"
    for completion in candidate_completion_stream(EFF_OPERATORS[witness.operator.name],
            witness.exotics,topology=witness.canonical_topology):
        decoded=loads_completion(dumps_completion(completion))
        if audit_amplitude_symmetrisation(decoded).status == "nonzero":
            break
    else:
        pytest.fail("Known copy-dependent mechanism was lost during construction")
    assert set(species_descriptors(decoded)) <= set(species_descriptors(witness))


@pytest.mark.parametrize("case",range(3))
def test_two_copy_witness_cannot_exclude_one_copy_candidate(case):
    witness=next(iter_completion_jsonl(FIXTURES/f"{case}.jsonl"))
    required=tuple(species_descriptors(witness).values())
    counts=Counter(required)
    repeated=next(key for key,count in counts.items() if count == 2)
    one_copy=list(required)
    one_copy.remove(repeated)
    with pytest.raises(ValueError,match="absent"):
        filter_species_registry({"candidate":{tuple(one_copy):None}}, {"candidate":1.0},
            mechanism_search=lambda *args:dict(required_species=required,verified=True),
            loop_search=lambda *args:None)


def test_candidate_derivative_search_retains_two_physical_numerators(monkeypatch):
    import importlib
    from neutrinomass.completions.operators import DERIV_EFF_OPERATORS
    from neutrinomass.completions.selection import generation_selection
    from neutrinomass.completions.routing import validate_routes
    module=importlib.import_module("neutrinomass.completions.completions")
    witness=next(iter_completion_jsonl(FIXTURES.parent/"physical-routes/D22-two-kinetic.jsonl"))
    operator=DERIV_EFF_OPERATORS["D22"]
    topologies=[x for x in module.get_topology_data(**operator.topology_type)
                if x["canonical_topology"] == witness.canonical_topology]
    monkeypatch.setattr(module,"get_topology_data",lambda **kwargs:topologies)
    with generation_selection(dict(producer="routed")):
        for completion in candidate_completion_stream(operator,witness.exotics):
            if sum(p.numerator_kind == "momentum" for p in completion.momentum_contributions) == 2:
                validate_routes(completion,derivative_order=2)
                assert audit_amplitude_symmetrisation(completion).status == "nonzero"
                break
        else:
            pytest.fail("Candidate-specific search lost a required two-numerator mechanism")
