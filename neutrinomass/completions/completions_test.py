#!/usr/bin/env python3

from neutrinomass.completions.completions import *
from neutrinomass.tensormethod.sm import L, Q, H, eb, ub, db
from neutrinomass.tensormethod.core import D

from neutrinomass.completions.operators import EFF_OPERATORS, DERIV_EFF_OPERATORS
from neutrinomass.completions.fingerprints import completion_digest

from sympy import Rational
import networkx as nx

from importlib import import_module
from pathlib import Path
from collections import Counter, defaultdict
from types import SimpleNamespace
import os
import subprocess
import sys

import pytest


def lnv_completions(op):
    return operator_completions(EFF_OPERATORS[op])


SIEVE = collect_completions(lnv_completions("1"))


def test_get_lorentz_epsilons():
    passes, epsilons = get_lorentz_epsilons((D(D(H, "11"), "00")("i0"), H("i1")))
    assert passes
    assert not epsilons

    passes, epsilons = get_lorentz_epsilons(
        (D(H, "11")("u0 d0 i0"), D(H, "11")("u1 d1 i1"))
    )
    assert passes
    assert epsilons == [eps("-u0 -u1"), eps("-d0 -d1")]

    passes, epsilons = get_lorentz_epsilons((Q("u0 c0 i0"), L("u1 i1")))
    assert passes
    assert epsilons == [eps("-u0 -u1")]

    passes, epsilons = get_lorentz_epsilons((Q("u0 c0 i0"), L.conj("d1 i1")))
    assert not passes
    assert not epsilons
    derivative_state = {"remaining": 1, "pending_fermions": [], "edges": []}
    assert not route_derivative_to_internal_fermion(
        (Q("u0 c0 i0"), L.conj("d1 i1")), derivative_state
    )
    assert derivative_state["remaining"] == 1

    passes, epsilons = get_lorentz_epsilons(
        (D(L, "01")("d0 i0"), D(H, "11")("u0 d1 i1"))
    )
    assert passes
    assert epsilons == [eps("-d0 -d1")]

    passes, epsilons = get_lorentz_epsilons(
        (D(Q, "01")("d0 c0 i0"), D(L, "01")("d1 i1"))
    )
    assert passes
    assert epsilons == [eps("-d0 -d1")]


def test_process_derivative_term():
    n = 10
    symbols = {
        "fermion": list(map(lambda i: "F" + str(i), range(n))),
        "boson": list(map(lambda i: "F" + str(i), range(n))),
    }
    exotic_dict = {}

    # Dirac fermion, deriv on fermion
    _, _, term = exotic_field_and_term(
        H("i1") * D(Q, "01")("d0 c0 i0"), symbols, exotic_dict
    )
    assert process_derivative_term(term)
    # Dirac fermion, deriv on scalar
    _, _, term = exotic_field_and_term(
        D(H, "11")("u0 d0 i1") * Q("u1 c0 i0") * eps("-u0 -u1"), symbols, exotic_dict
    )
    assert process_derivative_term(term)

    # Majorana fermion
    _, _, term = exotic_field_and_term(
        D(H, "11")("u0 d0 i1") * L("u1 i0") * eps("-u0 -u1"), symbols, exotic_dict
    )
    assert process_derivative_term(term)

    # Two derivatives, scalar case
    _, _, term = exotic_field_and_term(
        D(H, "11")("u0 d0 i1")
        * D(H, "11")("u1 d1 i0")
        * eps("-u0 -u1")
        * eps("-d0 -d1"),
        symbols,
        exotic_dict,
    )
    assert process_derivative_term(term)

    # Two derivatives, scalar case, four point
    _, _, term = exotic_field_and_term(
        D(H, "11")("u0 d0 i1")
        * D(H, "11")("u1 d1 i0")
        * H("i2")
        * eps("-u0 -u1")
        * eps("-d0 -d1"),
        symbols,
        exotic_dict,
    )
    assert process_derivative_term(term)

    # Two derivatives, fermion case
    _, _, term = exotic_field_and_term(
        D(L, "01")("d0 i1") * D(L, "01")("d1 i0") * eps("-d0 -d1"), symbols, exotic_dict
    )
    assert process_derivative_term(term)

    # Regular fermion case, for comparison above (and check on exotic_dict)
    _, _, term = exotic_field_and_term(
        L("u0 i1") * L("u1 i0") * eps("-u0 -u1"), symbols, exotic_dict
    )
    assert process_derivative_term(term)


def test_distinct_contraction_channels_get_distinct_exotic_labels():
    symbols = {"fermion": ["ψ", "χ"], "boson": ["φ", "η"]}
    field_dict = {}

    triplet, *_ = contract(
        (L("u0 i0"), H("i1")), symbols, [], field_dict
    )
    singlet, *_ = contract(
        (L("u1 i2"), H("i3")),
        symbols,
        [eps("-i2 -i3")],
        field_dict,
    )

    assert triplet.dynkin == "10002"
    assert singlet.dynkin == "10000"
    assert triplet.label != singlet.label


def test_construct_completion():
    data = partitions(EFF_OPERATORS["2"])[0]
    result = construct_completion(data["partition"], data["epsilons"], data["graph"])

    assert not isinstance(result, str)
    assert len(result) == 5


def test_completion_is_hashable():
    completion = next(operator_completions(EFF_OPERATORS["1"]))

    assert isinstance(hash(completion), int)


def test_exact_completion_mode_matches_raw_equivalence_classes():
    raw = list(operator_completions(EFF_OPERATORS["1"]))
    expected = []
    append_unique_completions(expected, raw)

    exact = exact_completions(EFF_OPERATORS["1"])

    assert len(exact) == len(expected)
    assert all(
        any(are_equivalent_completions(left, right) for right in exact)
        for left in expected
    )


def test_completions_dispatches_on_derivatives(monkeypatch):
    completions_module = import_module("neutrinomass.completions.completions")
    regular_result = object()
    derivative_result = object()
    monkeypatch.setattr(
        completions_module,
        "operator_completions",
        lambda operator, **kwargs: regular_result,
    )
    monkeypatch.setattr(
        completions_module,
        "deriv_operator_completions",
        lambda operator, **kwargs: derivative_result,
    )

    assert completions_module.completions(EFF_OPERATORS["1"]) is regular_result
    assert (
        completions_module.completions(DERIV_EFF_OPERATORS["D3"])
        is derivative_result
    )


def test_completions():
    # completions of dimension 5 and 7 operators
    o2 = collect_completions(lnv_completions("2"))
    o3a = collect_completions(lnv_completions("3a"))
    o3b = collect_completions(lnv_completions("3b"))
    o4a = collect_completions(lnv_completions("4a"))
    o4b = collect_completions(lnv_completions("4b"))
    o8 = collect_completions(lnv_completions("8"))

    # test on results of 1410.0689
    o2_comps = filter_completions(o2, SIEVE)
    o3a_comps = filter_completions(o3a, SIEVE)
    o3b_comps = filter_completions(o3b, SIEVE)
    o4a_comps = filter_completions(o4a, SIEVE)
    o4b_comps = filter_completions(o4b, SIEVE)
    o8_comps = filter_completions(o8, SIEVE)

    assert len(o2_comps) == 2
    assert len(o3a_comps) == 6
    assert len(o3b_comps) == 5
    assert not o4a_comps
    assert len(o4b_comps) == 3
    assert len(o8_comps) == 4


def test_collect_completions_preserves_nonadjacent_groups():
    completions = list(lnv_completions("1"))
    expected = defaultdict(list)
    for completion in completions:
        key = tuple(sorted(set(completion.exotic_info().values())))
        expected[key].append(completion)

    groups = list(expected.values())
    interleaved = [
        group[index]
        for index in range(max(map(len, groups)))
        for group in groups
        if index < len(group)
    ]
    collected = collect_completions(interleaved)

    assert {key: len(value) for key, value in collected.items()} == {
        key: len(value) for key, value in expected.items()
    }
    assert sum(map(len, collected.values())) == len(completions) == 8


def test_collect_completions_is_hash_seed_independent(tmp_path):
    script = """
from neutrinomass.completions import EFF_OPERATORS
from neutrinomass.completions.completions import collect_completions, operator_completions

completions = list(operator_completions(EFF_OPERATORS["1"]))
collected = collect_completions(completions)
print(len(completions), sum(map(len, collected.values())))
"""
    env = {**os.environ, "PYTHONHASHSEED": "1", "MPLCONFIGDIR": str(tmp_path)}
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=Path(__file__).resolve().parents[2],
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )

    assert result.stdout.strip() == "8 8"


def test_partition_order_is_hash_seed_independent(tmp_path):
    script = """
import hashlib
from neutrinomass.completions import EFF_OPERATORS
from neutrinomass.completions.completions import partitions

operator_partitions = partitions(EFF_OPERATORS["1"])
data = "\\n".join(repr(partition["partition"]) for partition in operator_partitions)
print(len(operator_partitions), hashlib.sha256(data.encode()).hexdigest())
"""
    outputs = []
    for seed in ("0", "1"):
        env = {
            **os.environ,
            "PYTHONHASHSEED": seed,
            "MPLCONFIGDIR": str(tmp_path / seed),
        }
        result = subprocess.run(
            [sys.executable, "-c", script],
            cwd=Path(__file__).resolve().parents[2],
            env=env,
            check=True,
            capture_output=True,
            text=True,
        )
        outputs.append(result.stdout.strip())

    assert outputs[0] == outputs[1]


def test_quick_remove_equivalent_partitions_preserves_first_occurrence():
    partitions_ = [("second",), ("first",), ("second",)]

    assert quick_remove_equivalent_partitions(partitions_) == [
        ("second",),
        ("first",),
    ]


def test_o9_completions():
    assert next(operator_completions(EFF_OPERATORS["9"]), None) is not None


def test_deriv_completions():
    # operators from the dimension-6 SMEFT
    from neutrinomass.tensormethod.sm import H, eb
    from neutrinomass.tensormethod.core import D

    Ophie_D1 = invariants(D(H, "11"), H.conj, eb.conj, eb)[0]
    Ophie_D2 = invariants(H, D(H.conj, "11"), eb.conj, eb)[0]
    Ophie_D3 = invariants(H, H.conj, D(eb.conj, "10"), eb)[0]
    Ophie_D4 = invariants(H, H.conj, eb.conj, D(eb, "01"))[0]

    eff_ophie_D1 = EffectiveOperator("OphieD1", Ophie_D1)
    eff_ophie_D2 = EffectiveOperator("OphieD2", Ophie_D2)
    eff_ophie_D3 = EffectiveOperator("OphieD3", Ophie_D3)
    eff_ophie_D4 = EffectiveOperator("OphieD4", Ophie_D4)

    models = {
        op.name: sorted(collect_completions(operator_completions(op)))
        for op in [eff_ophie_D1, eff_ophie_D2, eff_ophie_D3, eff_ophie_D4]
    }

    assert models["OphieD1"] == models["OphieD2"]
    assert models["OphieD1"] == models["OphieD3"]
    assert models["OphieD1"] == models["OphieD4"]


def test_d20_routes_derivative_to_internal_fermion(monkeypatch):
    completions_module = import_module("neutrinomass.completions.completions")
    topology_data = get_topology_data(5, 2)
    topology_11 = [
        data
        for data in topology_data
        if Path(data["partition_file"]).stem == "5s2f_11"
    ]
    monkeypatch.setattr(
        completions_module, "get_topology_data", lambda **kwargs: topology_11
    )

    completions = deriv_operator_completions(DERIV_EFF_OPERATORS["D20"])
    expected = tuple(
        sorted(
            (
                ("F", 0, 0, 3, ("3b", 0), ("y", Rational("1/2"))),
                ("S", 0, 0, 2, ("3b", 0), ("y", 0)),
                ("S", 0, 0, 3, ("3b", 0), ("y", Rational("3/2"))),
            )
        )
    )

    matching = [
        c for c in completions if tuple(sorted(c.exotic_info().values())) == expected
    ]
    forbidden_doublet_hhh = tuple(
        sorted(
            (
                ("F", 0, 0, 0, ("3b", 0), ("y", 1)),
                ("S", 0, 0, 0, ("3b", 0), ("y", 0)),
                ("S", 0, 0, 1, ("3b", 0), ("y", Rational("3/2"))),
            )
        )
    )

    assert matching
    assert not any(
        tuple(sorted(completion.exotic_info().values())) == forbidden_doublet_hhh
        for completion in completions
    )
    assert all(
        term.safe_simplify() != 0
        for completion in completions
        for term in completion.terms
    )
    assert all(len(c.derivative_edges) == 1 for c in matching)
    assert all(len(c.derivative_routes) == 1 for c in matching)
    assert all(c.topology == "5s2f_3" for c in completions)
    assert all(c.canonical_topology == "5s2f_11" for c in completions)

    julian = matching[0]
    assert {
        tuple(sorted(stringify_qns(field) for field in term.fields))
        for term in julian.terms
    } == {
        ("H", "H", "H", "S,00,3,-3/2,0"),
        ("F,00,3,-1/2,0", "S,00,3,3/2,0", "eb.conj"),
        ("H", "H.conj", "S,00,2,0,0"),
        ("F,00,3,1/2,0", "L", "S,00,2,0,0"),
    }
    assert all(
        is_singlet(term) and sum(field.mass_dim for field in term.fields) <= 4
        for term in julian.terms
    )
    assert julian.derivative_routes[0].differentiated_field == "L"


def test_derivative_routing_branches_and_is_order_independent(monkeypatch):
    completions_module = import_module("neutrinomass.completions.completions")
    topology_data = get_topology_data(5, 2)
    topology_19 = [
        data
        for data in topology_data
        if Path(data["partition_file"]).stem == "5s2f_19"
    ]
    monkeypatch.setattr(
        completions_module, "get_topology_data", lambda **kwargs: topology_19
    )

    original_partitions = completions_module.partitions

    def first_partition(operator, verbose=False):
        return original_partitions(operator, verbose=verbose)[:1]

    monkeypatch.setattr(completions_module, "partitions", first_partition)
    baseline = completions_module.momentum_routed_completions(
        DERIV_EFF_OPERATORS["D20"]
    )

    assert len(baseline) == 2
    assert all(len(completion.derivative_routes) == 1 for completion in baseline)
    assert all(len(completion.derivative_edges) == 1 for completion in baseline)
    assert len({completion.derivative_edges[0] for completion in baseline}) == 2

    original_candidates = completions_module.derivative_route_candidates
    monkeypatch.setattr(
        completions_module,
        "derivative_route_candidates",
        lambda fields: tuple(reversed(original_candidates(fields))),
    )
    reversed_candidates = completions_module.momentum_routed_completions(
        DERIV_EFF_OPERATORS["D20"]
    )

    def reverse_children(partition):
        if isinstance(partition, Leaf):
            return partition
        return tuple(reverse_children(branch) for branch in reversed(partition))

    def rerooted_partition(operator, verbose=False):
        partition = original_partitions(operator, verbose=verbose)[0]
        roots = canonical_rooted_partitions(
            partition["partition"], partition["graph"]
        )
        partition = dict(partition)
        partition["partition"] = reverse_children(roots[-1])
        return [partition]

    monkeypatch.setattr(completions_module, "partitions", rerooted_partition)
    rerooted = completions_module.momentum_routed_completions(
        DERIV_EFF_OPERATORS["D20"]
    )

    for comparison in (reversed_candidates, rerooted):
        assert all(
            any(are_equivalent_completions(left, right) for right in comparison)
            for left in baseline
        )
        assert all(
            any(are_equivalent_completions(right, left) for left in baseline)
            for right in comparison
        )


def test_derivative_rerooting_enumerates_every_internal_vertex():
    graph = nx.Graph([(0, 1), (1, 2), (0, 3), (2, 4)])
    partition = (Leaf("S", 3), Leaf("F", 4))

    rooted = canonical_rooted_partitions(partition, graph)

    assert len(rooted) == 3
    assert all(
        sorted(leaf.node for leaf in partition_leaves(candidate)) == [3, 4]
        for candidate in rooted
    )


def test_derivative_route_choice_count_uses_fermion_parity():
    scalar = SimpleNamespace(is_fermion=False)
    fermion = SimpleNamespace(is_fermion=True)

    no_internal_fermion = (Leaf(fermion, 0), Leaf(fermion, 1), Leaf(scalar, 2))
    one_internal_fermion = (
        (Leaf(fermion, 0), Leaf(scalar, 1)),
        Leaf(fermion, 2),
        Leaf(scalar, 3),
    )
    two_internal_fermions = (
        (Leaf(fermion, 0), Leaf(scalar, 1)),
        (Leaf(fermion, 2), Leaf(scalar, 3)),
        Leaf(scalar, 4),
    )

    assert derivative_route_choice_count(no_internal_fermion) == 0
    assert derivative_route_choice_count(one_internal_fermion) == 1
    assert derivative_route_choice_count(two_internal_fermions) == 2


def test_derivative_routing_support_is_explicit():
    supported = {
        name
        for name, operator in DERIV_EFF_OPERATORS.items()
        if operator_strip_derivs(operator.operator)["n_derivs"] == 1
        and (
            unique_lorentz_completion_operator(operator) is not None
            or name in PROJECTED_LORENTZ_OPERATORS
        )
    }
    skipped = {
        name
        for name, operator in DERIV_EFF_OPERATORS.items()
        if operator_strip_derivs(operator.operator)["n_derivs"] == 1
        and unique_lorentz_completion_operator(operator) is None
        and name not in PROJECTED_LORENTZ_OPERATORS
    }

    assert supported == {
        "D3",
        "D5a",
        "D5b",
        "D5c",
        "D5d",
        "D10a",
        "D10b",
        "D10c",
        "D20",
    }
    assert skipped == {
        "D6a",
        "D6b",
        "D8a",
        "D8b",
        "D8c",
        "D8d",
        "D8e",
        "D8f",
        "D8g",
        "D8h",
        "D8i",
        "D9a",
        "D9b",
        "D12a",
        "D12b",
        "D14a",
        "D14b",
        "D14c",
        "D16a",
        "D16b",
        "D16c",
        "D17",
    }
    assert momentum_routed_completions(DERIV_EFF_OPERATORS["D6a"]) == []


def test_d6a_projected_routing_has_a_positive_uv_model(monkeypatch):
    completions_module = import_module("neutrinomass.completions.completions")
    original_partitions = completions_module.partitions
    monkeypatch.setattr(
        completions_module, "PROJECTED_LORENTZ_OPERATORS", frozenset({"D6a"})
    )

    def projected_partition(operator, verbose=False):
        return original_partitions(operator, verbose=verbose)[144:145]

    monkeypatch.setattr(completions_module, "partitions", projected_partition)
    completions = completions_module.momentum_routed_completions(
        DERIV_EFF_OPERATORS["D6a"]
    )

    assert len(completions) == 4
    assert all(
        completion.lorentz_projection.coordinates == ("1", "-1")
        for completion in completions
    )
    assert all(
        completion.lorentz_projection.basis_labels
        == (
            "u:0-1,2-3|d:0-1",
            "u:0-2,1-3|d:0-1",
        )
        for completion in completions
    )
    assert all(
        completion_quantum_number_strings(completion)
        == {
            "F,00,1,1/2,0",
            "F,00,1,3/2,0",
            "S,00,2,1,0",
        }
        for completion in completions
    )
    assert all(
        is_singlet(term) and term.safe_simplify() != 0
        for completion in completions
        for term in completion.terms
    )


def test_d6a_ibp_eom_cut_weights_cover_equivalent_and_zero_projections():
    operator = DERIV_EFF_OPERATORS["D6a"]
    fields, epsilons, _ = operator_strip_derivs(operator.operator).values()
    stripped = construct_operator(fields, epsilons)
    derivative = next(
        field for field in operator.operator.indexed_fields if field.derivs
    )
    target_gauge = tuple(map(str, derivative.gauge_indices))
    target = next(
        field
        for field in stripped.indexed_fields
        if field.field == derivative.strip_derivs()
        and tuple(map(str, field.gauge_indices)) == target_gauge
    )
    other = next(
        field
        for field in stripped.indexed_fields
        if field.label == target.label
        and field.dynkin == target.dynkin
        and field.charges == target.charges
        and field is not target
    )
    spectators = [
        field for field in stripped.indexed_fields if field not in (target, other)
    ]
    route = DerivativeRoute((6, 7), "F", "10", "L", "10")
    graph = nx.Graph(
        [(6, 7), (6, 0), (6, 1), (6, 2), (7, 3), (7, 4), (7, 5)]
    )

    def weight(left_higgs, right_higgs):
        ordered = [left_higgs, *spectators[:2], right_higgs, *spectators[2:]]
        partition = tuple(Leaf(field, node) for node, field in enumerate(ordered))
        return routed_ibp_weight(operator.operator, partition, graph, route)

    assert weight(target, other) == 1
    assert weight(other, target) == -1

    ordered = [target, other, *spectators]
    partition = tuple(Leaf(field, node) for node, field in enumerate(ordered))
    assert routed_ibp_weight(operator.operator, partition, graph, route) == 0


def test_operator_strip_derivs_preserves_field_statistics():
    stripped = operator_strip_derivs(DERIV_EFF_OPERATORS["D20"].operator)
    higgs_fields = [
        field for field, _ in stripped["fields"] if field.label == H.label
    ]

    assert len(higgs_fields) == 5
    assert all(field.comm == H.comm for field in higgs_fields)
    assert all(field.nf == H.nf for field in higgs_fields)


def test_d3_full_census_baseline_is_stable():
    completions = deriv_operator_completions(DERIV_EFF_OPERATORS["D3"])
    unique = []
    append_unique_completions(unique, completions)

    assert len(completions) == 26
    assert sum(bool(item.derivative_routes) for item in completions) == 2
    assert Counter(
        (item.topology, item.canonical_topology) for item in completions
    ) == {
        ("3s2f_3", "3s2f_3"): 18,
        ("3s2f_4", "3s2f_4"): 8,
    }
    assert completion_digest(completions) == (
        "dd82599a00ecd0ed202fa9c677e66ae06801828fa7b2b723f275e3515f5e1d4b"
    )

    assert len(unique) == 7
    assert sum(bool(item.derivative_routes) for item in unique) == 2
    assert completion_digest(unique) == (
        "db7f3b19d4a595b62bf008a98829030e9ae44c199764c70217bba39fb4db7ebf"
    )

    canonical = exact_completions(DERIV_EFF_OPERATORS["D3"])
    assert len(canonical) == len(unique)
    assert completion_digest(canonical) == completion_digest(unique)


ROUTED_MODEL_CONTROLS = {
    "D3": ("3s2f_3", {"F,00,1,1/2,0", "F,00,2,0,0"}),
    "D5a": (
        "2s4f_4",
        {"F,00,2,0,0", "F,00,2,1,0", "S,00,2,1,0"},
    ),
    "D5b": (
        "2s4f_4",
        {"F,00,0,0,0", "F,00,0,1,0", "S,00,0,1,0"},
    ),
    "D5c": (
        "2s4f_4",
        {"F,00,2,0,0", "F,00,2,1,0", "S,00,2,1,0"},
    ),
    "D5d": (
        "2s4f_4",
        {"F,00,0,0,0", "F,00,2,1,0", "S,00,2,1,0"},
    ),
    "D10a": (
        "2s4f_4",
        {"F,00,2,0,0", "F,10,1,1/6,1", "S,10,1,1/6,1"},
    ),
    "D10b": (
        "2s4f_4",
        {"F,00,2,0,0", "F,10,1,1/6,1", "S,10,1,1/6,1"},
    ),
    "D10c": (
        "2s4f_4",
        {"F,00,2,1,0", "F,10,1,7/6,1", "S,10,1,1/6,1"},
    ),
    "D20": (
        "5s2f_10",
        {"F,00,2,0,0", "F,00,3,1/2,0", "S,00,3,3/2,0"},
    ),
}


def completion_quantum_number_strings(completion):
    quantum_numbers = set()
    for lorentz, colour_up, colour_down, isospin, (_, baryon), (_, hypercharge) in (
        completion.exotic_info().values()
    ):
        quantum_numbers.add(
            f"{lorentz},{colour_up}{colour_down},{isospin},{hypercharge},{baryon}"
        )
    return quantum_numbers


@pytest.mark.parametrize("operator_name", ROUTED_MODEL_CONTROLS)
def test_supported_routing_has_physics_controls(operator_name, monkeypatch):
    completions_module = import_module("neutrinomass.completions.completions")
    topology, positive_model = ROUTED_MODEL_CONTROLS[operator_name]
    operator = DERIV_EFF_OPERATORS[operator_name]
    topology_data = [
        data
        for data in get_topology_data(**operator.topology_type)
        if data["canonical_topology"] == topology
    ]
    monkeypatch.setattr(
        completions_module, "get_topology_data", lambda **kwargs: topology_data
    )

    completions = completions_module.momentum_routed_completions(operator)
    matching = [
        completion
        for completion in completions
        if completion_quantum_number_strings(completion) == positive_model
    ]

    assert matching
    assert all(len(completion.derivative_routes) == 1 for completion in matching)
    assert all(
        is_singlet(term)
        and term.safe_simplify() != 0
        and sum(field.mass_dim for field in term.fields) <= 4
        for completion in matching
        for term in completion.terms
    )
    for completion in matching:
        route = completion.derivative_routes[0]
        particle = completion.graph.edges[route.edge]["particle"]
        assert base_exotic_label(particle) == base_exotic_label(
            route.numerator_field
        )

    term = matching[0].terms[0]
    field = term.indexed_fields[0]
    shifted_charges = dict(field.charges)
    shifted_charges["y"] += 1
    non_singlet_field = IndexedField(
        field.label,
        " ".join(str(index) for index in field.indices),
        charges=shifted_charges,
        is_conj=field.is_conj,
        comm=field.comm,
    )
    non_singlet = Operator(
        *(non_singlet_field if tensor is field else tensor for tensor in term.tensors)
    )
    assert not is_singlet(non_singlet)


def test_derivs_nlo_completions():
    # Derivative operator examples:
    from neutrinomass.tensormethod.sm import L, H
    from neutrinomass.tensormethod.core import D

    o1box = EffectiveOperator(
        "O1box",
        L("u0 i0")
        * L("u1 i1")
        * D(H, "11")("u2 d0 i2")
        * D(H, "11")("u3 d1 i3")
        * eps("-i0 -i2")
        * eps("-i1 -i3"),
    )

    comps = collect_completions(operator_completions(o1box))

    assert len(comps) == 3
    assert not filter_completions(comps, SIEVE)


def test_ophibox_ophiD():
    """Example from section 2.2 in the paper."""

    from neutrinomass.tensormethod.sm import H
    from neutrinomass.tensormethod.core import D

    ohhdd1 = EffectiveOperator(
        "OHHDD1",
        H.conj("i0")
        * H.conj("i1")
        * D(H, "11")("u0 d0 i2")
        * D(H, "11")("u1 d1 i3")
        * eps("-i0 -i2")
        * eps("-i1 -i3"),
    )

    ohhdd2 = EffectiveOperator(
        "OHHDD2",
        H.conj("i0")
        * H("i1")
        * D(H.conj, "11")("u0 d0 i2")
        * D(H, "11")("u1 d1 i3")
        * eps("-i0 -i1")
        * eps("-i2 -i3"),
    )

    ohhdd3 = EffectiveOperator(
        "OHHDD3",
        H.conj("i0")
        * H("i1")
        * D(H.conj, "11")("u0 d0 i2")
        * D(H, "11")("u1 d1 i3")
        * eps("-i0 -i2")
        * eps("-i1 -i3"),
    )

    ohhdd4 = EffectiveOperator(
        "OHHDD4",
        H.conj("i0")
        * H("i1")
        * D(H.conj, "11")("u0 d0 i2")
        * D(H, "11")("u1 d1 i3")
        * eps("-i0 -i3")
        * eps("-i1 -i2"),
    )

    comps = {}
    for op in [ohhdd1, ohhdd2, ohhdd3, ohhdd4]:
        comps[op.name] = list(collect_completions(operator_completions(op)))

    assert len(comps["OHHDD1"]) == 1
    assert len(comps["OHHDD2"]) == 1
    assert len(comps["OHHDD3"]) == 1
    assert len(comps["OHHDD4"]) == 1

    assert comps["OHHDD1"][0] == (("S", 0, 0, 2, ("3b", 0), ("y", 1)),)
    assert comps["OHHDD2"][0] == (("S", 0, 0, 0, ("3b", 0), ("y", 0)),)
    assert comps["OHHDD3"][0] == (("S", 0, 0, 2, ("3b", 0), ("y", 0)),)
    assert comps["OHHDD4"][0] == (("S", 0, 0, 2, ("3b", 0), ("y", 0)),)


def test_1204_5986_completions():
    """Paper 1204.5986 lists UV completions of some derivative operators. Check
    output of program against these results.

    """

    from neutrinomass.tensormethod.sm import eb, H, L
    from neutrinomass.tensormethod.core import D

    # fields
    k = ("S", 0, 0, 0, ("3b", 0), ("y", 2))
    xi1 = ("S", 0, 0, 2, ("3b", 0), ("y", 1))
    sigma = ("F", 0, 0, 2, ("3b", 0), ("y", 0))
    ltilde = ("F", 0, 0, 1, ("3b", 0), ("y", Rational("1/2")))
    z = ("S", 0, 0, 1, ("3b", 0), ("y", Rational("3/2")))
    nur = ("F", 0, 0, 0, ("3b", 0), ("y", 0))

    # their notation
    phi_0_2 = ("S", 0, 0, 0, ("3b", 0), ("y", 2))
    phi_12_32 = ("S", 0, 0, 1, ("3b", 0), ("y", Rational("3/2")))
    phi_1_1 = ("S", 0, 0, 2, ("3b", 0), ("y", 1))
    psi_1_0 = ("F", 0, 0, 2, ("3b", 0), ("y", 0))
    psi_12_12 = ("F", 0, 0, 1, ("3b", 0), ("y", Rational("1/2")))
    psi_0_0 = ("F", 0, 0, 0, ("3b", 0), ("y", 0))

    o9 = EffectiveOperator(
        "O9",
        D(H, "11")("u0 d0 i0")
        * D(H, "11")("u1 d1 i1")
        * H("i2")
        * H("i3")
        * eb.conj("d2 g0")
        * eb.conj("d3 g1")
        * eps("-i0 -i2")
        * eps("-i1 -i3"),
    )

    # models from table 5 in 1204.5986 (the ones implying vanishing vertices
    # have been left out)
    o9_models_from_paper = {
        # first topology
        frozenset([phi_0_2, phi_12_32, phi_1_1]),
        frozenset([psi_12_12, phi_12_32, phi_1_1]),
        frozenset([psi_12_12, psi_1_0, phi_1_1]),
        frozenset([psi_12_12, psi_0_0]),
        frozenset([psi_12_12, psi_1_0]),
        frozenset([phi_1_1, psi_1_0]),
        # second topology
        frozenset([phi_0_2, phi_1_1]),
        frozenset([psi_12_12, phi_1_1]),
    }

    o7_models_from_paper = {
        frozenset([psi_1_0, phi_1_1]),
        frozenset([psi_12_12, phi_1_1]),
        frozenset([psi_12_12, psi_0_0]),
        frozenset([psi_12_12, psi_1_0]),
    }

    o7 = EffectiveOperator(
        "O7",
        eb.conj("d0 g0")
        * L("u0 i0")
        * D(H, "11")("u1 d1 i1")
        * H("i2")
        * H("i3")
        * eps("-i0 -i2")
        * eps("-i1 -i3"),
    )

    o7_comps = collect_completions(operator_completions(o7))
    o9_comps = collect_completions(operator_completions(o9))
    out = []
    for comps, models in [
        (o7_comps, o7_models_from_paper),
        (o9_comps, o9_models_from_paper),
    ]:
        for k, v in comps.items():
            # we have, they don't have
            if frozenset(k) not in models:
                print("They don't have:")
                print((k, v))
                out.append((k, v))

        for model in models:
            # they have, we don't have
            if model not in set(map(frozenset, comps.keys())):
                print("We don't have:")
                print(model)

    assert not out
    # return out


def test_ibp():
    from itertools import combinations
    from collections import defaultdict

    od12_1 = EffectiveOperator(
        "OD12_1",
        D(L, "01")("d3 i0 g0")
        * Q("u1 c1 i1 g1")
        * eb.conj("d0 g2")
        * db("u2 -c2")
        * H("i2")
        * H("i3")
        * eps("-i0 -i2")
        * eps("-i1 -i3"),
    )
    od12_2 = EffectiveOperator(
        "OD12_2",
        L("u0 i0 g0")
        * D(Q, "01")("d3 c1 i1 g1")
        * eb.conj("d0 g2")
        * db("u2 -c2")
        * H("i2")
        * H("i3")
        * eps("-i0 -i2")
        * eps("-i1 -i3"),
    )
    od12_3 = EffectiveOperator(
        "OD12_3",
        L("u0 i0 g0")
        * Q("u1 c1 i1 g1")
        * D(eb.conj, "10")("u5 g2")
        * db("u2 -c2")
        * H("i2")
        * H("i3")
        * eps("-i0 -i2")
        * eps("-i1 -i3"),
    )
    od12_4 = EffectiveOperator(
        "OD12_4",
        L("u0 i0 g0")
        * Q("u1 c1 i1 g1")
        * eb.conj("d0 g2")
        * D(db, "01")("d5 -c2")
        * H("i2")
        * H("i3")
        * eps("-i0 -i2")
        * eps("-i1 -i3"),
    )
    od12_5 = EffectiveOperator(
        "OD12_5",
        L("u0 i0 g0")
        * Q("u1 c1 i1 g1")
        * eb.conj("d0 g2")
        * db("u2 -c2")
        * D(H, "11")("u3 d1 i2")
        * H("i3")
        * eps("-i0 -i2")
        * eps("-i1 -i3"),
    )
    od12_6 = EffectiveOperator(
        "OD12_6",
        L("u0 i0 g0")
        * Q("u1 c1 i1 g1")
        * eb.conj("d0 g2")
        * db("u2 -c2")
        * H("i2")
        * D(H, "11")("u3 d1 i3")
        * eps("-i0 -i2")
        * eps("-i1 -i3"),
    )

    models = {
        op.name: sorted(collect_completions(operator_completions(op)))
        for op in [od12_1, od12_2, od12_3, od12_4, od12_5, od12_6]
    }

    model_names = list(map(lambda i: "OD12_" + str(i), range(1, 7)))
    model_dict = defaultdict(set)
    for a, b in combinations(model_names, 2):
        for model in models[a]:
            if model in models[b]:
                model_dict[model].add(a)
                model_dict[model].add(b)

    represented_operators = set().union(*model_dict.values())
    assert represented_operators == set(model_names)
    assert all(len(names) >= 2 for names in model_dict.values())


def test_symmetries():
    o2_models = collect_models(lnv_completions("2"))
    o3a_models = collect_models(lnv_completions("3a"))
    o3b_models = collect_models(lnv_completions("3b"))
    o8_models = collect_models(lnv_completions("8"))

    for m in [*o2_models, *o3a_models, *o3b_models, *o8_models]:
        for c in m.completions:
            assert c.lagrangian.num_u1_symmetries() == 2


def test_compare_terms():
    phi1 = IndexedField("φ", "c0 i0 i1", charges={"y": -Rational(1, 3), "3b": 1})
    eta1 = IndexedField("η", "-c1 i2", charges={"y": -Rational(1, 6), "3b": -1})
    omega1 = IndexedField("ω", "i3", charges={"y": Rational(1, 2), "3b": 0})
    terms_1 = [
        phi1.conj
        * L("u0 i2 g378_")
        * Q("u1 c0 i3 g381_")
        * eps("-u0 -u1")
        * eps("-i0 -i2")
        * eps("-i1 -i3"),
        eta1.conj
        * L("u0 i1 g379_")
        * db("u1 -c1 g383_")
        * eps("-u0 -u1")
        * eps("-i2 -i1"),
        omega1.conj
        * L("u0 i1 g380_")
        * eb("u1 g385_")
        * eps("-u0 -u1")
        * eps("-i3 -i1"),
        phi1 * eta1 * omega1 * eps("-i0 -i3") * eps("-i1 -i2") * delta("c1 -c0"),
    ]

    omega2 = IndexedField("ω", "c0 i0 i1", charges={"y": -Rational(1, 3), "3b": 1})
    phi2 = IndexedField("φ", "-c1 i2", charges={"y": -Rational(1, 6), "3b": -1})
    eta2 = IndexedField("η", "i3", charges={"y": Rational(1, 2), "3b": 0})
    terms_2 = [
        phi2.conj
        * db("u0 -c1 g383_")
        * L("u1 i1 g379_")
        * eps("-u0 -u1")
        * eps("-i2 -i1"),
        eta2.conj * L("u0 i1 g378_") * eb("u1 g385_") * eps("-u0 -u1") * eps("-i3 -i1"),
        omega2.conj
        * Q("u0 c0 i2 g381_")
        * L("u1 i3 g380_")
        * eps("-u0 -u1")
        * eps("-i0 -i3")
        * eps("-i1 -i2"),
        phi2 * eta2 * omega2 * eps("-i1 -i2") * eps("-i3 -i0") * delta("c1 -c0"),
    ]

    remapping = {"φ": "ω", "η": "φ", "ω": "η"}
    assert check_remapping_on_terms(terms_1, terms_2, remapping) == remapping


def test_remove_duplicate_completions():
    comps = list(operator_completions(EFF_OPERATORS["3b"]))
    comps_len = len(comps)
    comps = clean_completions(comps)
    assert comps_len > len(comps)
