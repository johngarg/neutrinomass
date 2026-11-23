import importlib
from pathlib import Path
import pytest

from neutrinomass.completions.completions import (
    HistoricalDerivativeBasis, UnreducedSecondDerivativeBasis, operator_completions,
    deriv_operator_completion_stream, partition_leaves,
)
from neutrinomass.completions.operators import EFF_OPERATORS, DERIV_EFF_OPERATORS
from neutrinomass.completions.physical_routes import (
    local_physical_expansions, sewn_fermion_operator, project_second_derivative_tree,
)
from neutrinomass.completions.routing import validate_routes
from neutrinomass.tensormethod.lorentz import lorentz_contraction
from neutrinomass.database.serialization import iter_completion_jsonl,loads_completion,dumps_completion


def test_sewn_mass_tree_has_closed_lorentz_contraction():
    for completion in operator_completions(EFF_OPERATORS["1"]):
        assert lorentz_contraction(sewn_fermion_operator(completion)).coefficient


def test_local_second_order_expansions_have_exact_budget_and_physical_endpoints(monkeypatch):
    module = importlib.import_module("neutrinomass.completions.completions")
    operator = DERIV_EFF_OPERATORS["D22"]
    topologies = module.get_topology_data(**operator.topology_type)[:1]
    monkeypatch.setattr(module,"get_topology_data",lambda **kwargs:topologies)
    basis = UnreducedSecondDerivativeBasis.from_operator(operator)
    local = [c for c in deriv_operator_completion_stream(operator)
             if any(leaf.field.derivs for leaf in partition_leaves(c.partition))]
    expanded = [candidate for c in local for candidate in local_physical_expansions(
        c,2,second_derivative_basis=basis)]
    assert expanded
    for completion in expanded:
        assert completion.momentum_contributions
        validate_routes(completion,derivative_order=2)
        assert any(x != "0" for x in completion.lorentz_projection.coordinates)


def test_two_kinetic_lines_match_explicit_weyl_chain_after_serialization():
    from neutrinomass.tensormethod.core import Operator,eps
    from neutrinomass.completions.amplitudes import audit_amplitude_symmetrisation
    fixture=Path(__file__).resolve().parents[2]/"tests/fixtures/physical-routes/D22-two-kinetic.jsonl"
    completion=loads_completion(dumps_completion(next(iter_completion_jsonl(fixture))))
    basis=UnreducedSecondDerivativeBasis.from_operator(DERIV_EFF_OPERATORS["D22"])
    sewn=sewn_fermion_operator(completion,basis.momentum_fields)
    external=[leaf.field.strip_derivs_with_indices() for leaf in partition_leaves(completion.partition)]
    fermions=[f for f in external if f.is_fermion]
    momenta=[field(f"utest{i} dtest{i}") for i,field in enumerate(basis.momentum_fields)]
    # The saved chain has two dotted external spinors and two kinetic
    # numerators connected through an undotted mass contraction.
    manual=Operator(*external,*momenta,
        eps(f"{-fermions[0].indices_by_type['Dotted'][0]} -dtest0"),
        eps(f"{-fermions[1].indices_by_type['Dotted'][0]} -dtest1"),
        eps("-utest0 -utest1"))
    assert lorentz_contraction(sewn,normalise=True) == lorentz_contraction(manual,normalise=True)
    assert project_second_derivative_tree(basis,completion) == completion.lorentz_projection
    assert audit_amplitude_symmetrisation(completion).status == "nonzero"
    validate_routes(completion,derivative_order=2)
    completion.momentum_contributions=()
    with pytest.raises(ValueError):
        validate_routes(completion,derivative_order=2)


def test_two_kinetic_branch_is_enumerated_without_a_forced_vector_current(monkeypatch):
    from neutrinomass.completions.completions import momentum_routed_completion_stream
    module=importlib.import_module("neutrinomass.completions.completions")
    topologies=[x for x in module.get_topology_data(**DERIV_EFF_OPERATORS["D22"].topology_type)
                if x["canonical_topology"] == "4s2f_6"]
    monkeypatch.setattr(module,"get_topology_data",lambda **kwargs:topologies)
    for completion in momentum_routed_completion_stream(DERIV_EFF_OPERATORS["D22"]):
        if sum(p.numerator_kind == "momentum" for p in completion.momentum_contributions) == 2:
            validate_routes(completion,derivative_order=2)
            break
    else:
        pytest.fail("Two kinetic lines were lost by the vertex-driven routing heuristic")
