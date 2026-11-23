"""Physical endpoint regressions for the October 2026 routing defect."""
from functools import reduce

import networkx as nx
import pytest

from neutrinomass.completions.completions import (
    derivative_route_alternatives, process_derivative_term, record_vertex_ports,
)
from neutrinomass.completions.core import (
    MajoranaFermion, VectorLikeDiracFermion, RealScalar,
)


@pytest.mark.parametrize("vectorlike", [False, True])
@pytest.mark.parametrize("reverse", [False, True])
def test_selected_numerator_is_endpoint_changed_by_processing(vectorlike, reverse):
    cls = VectorLikeDiracFermion if vectorlike else MajoranaFermion
    charges = {"y": 1} if vectorlike else {"y": 0}
    left = cls("A", "u0", charges=charges)
    right = cls("B", "u1", charges=charges).conj
    fields = (left, right, RealScalar("S", ""))
    if reverse:
        fields = fields[::-1]
    alternatives = derivative_route_alternatives(fields)
    assert len(alternatives) == 2
    for routed, epsilons, (numerator, helper) in alternatives:
        term = process_derivative_term(reduce(lambda a, b: a * b, (*routed, *epsilons)))
        assert not isinstance(term, str)
        # Follow explicit partner operations, including the conjugate gauge
        # representation, rather than accepting a matching edge/particle label.
        changed = numerator.dirac_partner().conj if vectorlike else numerator.conj
        assert any(f.field == changed.field for f in term.indexed_fields)
        assert any((f.label, f.dynkin, f.charges, f.is_conj) ==
                   (helper.label, helper.dynkin, helper.charges, helper.is_conj)
                   for f in term.indexed_fields)
        assert numerator is not helper


@pytest.mark.parametrize("vectorlike", [False, True])
def test_one_heavy_one_light_uses_light_helper(vectorlike):
    from neutrinomass.tensormethod.core import IndexedField
    cls = VectorLikeDiracFermion if vectorlike else MajoranaFermion
    heavy = cls("N", "u0", charges={"y": 1 if vectorlike else 0})
    light = IndexedField("light", "d0", charges={"y": -heavy.y}, comm=1)
    routes = derivative_route_alternatives((heavy, light, RealScalar("S", "")))
    assert len(routes) == 1
    fields, epsilons, (numerator, helper) = routes[0]
    assert numerator is heavy and helper is light
    assert not isinstance(process_derivative_term(reduce(lambda a, b: a*b, (*fields, *epsilons))), str)


def test_vectorlike_pair_requires_weyl_partner_and_gauge_conjugation():
    from neutrinomass.completions.routing import species_states, classify_fermion_pair, InvalidMomentumGraph
    field = VectorLikeDiracFermion("F", "u0 c0", charges={"y": 1})
    states = species_states({field})
    assert classify_fermion_pair(field, field.conj, states) == "momentum"
    assert classify_fermion_pair(field, field.dirac_partner(), states) == "mass"
    with pytest.raises(InvalidMomentumGraph):
        classify_fermion_pair(field, field, states)
    with pytest.raises(InvalidMomentumGraph):
        classify_fermion_pair(field, field.dirac_partner().conj, states)


def saved_pair():
    from pathlib import Path
    from neutrinomass.database.serialization import loads_completion
    path = Path(__file__).parents[2] / "tests/fixtures/scalar-n-routing.jsonl"
    return [loads_completion(line) for line in path.read_text().splitlines()]


def test_wrong_edge_with_matching_particle_rejected_by_consumers():
    from neutrinomass.completions.routing import InvalidMomentumGraph
    from neutrinomass.completions.amplitudes import derivative_placement_polynomial
    from neutrinomass.completions.fingerprints import completion_fingerprint
    from neutrinomass.completions.completions import are_equivalent_completions
    from neutrinomass.database.serialization import dumps_completion
    from neutrinomass.database.export import export_completion
    corrupted = saved_pair()[0]
    route = corrupted.derivative_routes[0]
    assert corrupted.graph.edges[route.edge]["particle"] == route.numerator_field
    for consume in [derivative_placement_polynomial, completion_fingerprint,
                    dumps_completion, export_completion,
                    lambda c: are_equivalent_completions(c, c)]:
        with pytest.raises(InvalidMomentumGraph, match="disagree"):
            consume(corrupted)


def test_repeated_singlet_ports_follow_the_original_factors():
    from neutrinomass.completions.routing import InvalidMomentumGraph
    from neutrinomass.tensormethod.core import eps

    left = MajoranaFermion("N", "u0")
    right = MajoranaFermion("N", "u1")
    term = right * left * eps("-u1 -u0")
    graph = nx.path_graph(3)
    edges = ((0, 1), (1, 2))
    record_vertex_ports(graph, 1, 0, term, (left, right), edges)
    for port, source, edge in zip(graph.nodes[1]["heavy_ports"], (left, right), edges):
        assert term.tensors[port["tensor"]] is source
        assert tuple(port["edge"]) == edge

    # With both original occurrences lost, the singlets have identical
    # quantum numbers and no spectator gauge indices to distinguish them.
    detached = (MajoranaFermion("N", "u2"), MajoranaFermion("N", "u3"))
    with pytest.raises(InvalidMomentumGraph, match="Unresolved equal heavy occurrences"):
        record_vertex_ports(graph, 1, 0, term, detached, edges)


def test_repeated_species_ports_and_relabelling_survive_serialization():
    import json
    from neutrinomass.database.serialization import completion_to_record, completion_from_record
    from neutrinomass.completions.routing import migrate_legacy_one_derivative, validate_routes, InvalidMomentumGraph
    record = completion_to_record(migrate_legacy_one_derivative(saved_pair()[0]))
    # Identify two Majorana species while retaining their distinct concrete
    # occurrence ports. Their adjacent equal species cannot be sorted by label.
    payload = json.dumps(record, ensure_ascii=False).replace('"f"', '"ψ"')
    repeated = completion_from_record(json.loads(payload))
    validate_routes(repeated, derivative_order=1)
    relabelled = json.loads(payload.replace('ψ', 'renamed'))
    validate_routes(completion_from_record(relabelled), derivative_order=1)
    for _, attrs in repeated.graph.nodes(data=True):
        attrs.pop("heavy_ports", None)
    with pytest.raises(InvalidMomentumGraph, match="repeated-species"):
        validate_routes(repeated)


def test_denominator_order_is_separate_from_endpoint_chirality():
    from neutrinomass.completions.routing import migrate_legacy_one_derivative, validate_routes, InvalidMomentumGraph
    completion = migrate_legacy_one_derivative(saved_pair()[1])
    completion.momentum_contributions = (completion.momentum_contributions[0]._replace(denominator_order=1),)
    validate_routes(completion, derivative_order=3)
    with pytest.raises(InvalidMomentumGraph, match="budget"):
        validate_routes(completion, derivative_order=1)


def test_saved_zero_and_nonzero_records_and_roundtrip():
    from pathlib import Path
    from neutrinomass.database.serialization import loads_completion, dumps_completion
    from neutrinomass.completions.routing import (
        inspect_endpoints, validate_routes, migrate_legacy_one_derivative, InvalidMomentumGraph,
    )
    from neutrinomass.completions.completions import partition_leaves
    path = Path(__file__).parents[2] / "tests/fixtures/scalar-n-routing.jsonl"
    raw = [loads_completion(line) for line in path.read_text().splitlines()]
    assert inspect_endpoints(raw[0]).momentum_edges == ((17496, 157464),)
    with pytest.raises(InvalidMomentumGraph, match="disagree"):
        validate_routes(raw[0], derivative_order=1)
    for original, expected in zip(raw, [0, 1]):
        repaired = migrate_legacy_one_derivative(original)
        validate_routes(repaired, derivative_order=1)
        restored = loads_completion(dumps_completion(repaired))
        validate_routes(restored, derivative_order=1)
        leaves = {leaf.node: leaf.field for leaf in partition_leaves(restored.partition)}
        # Restricted constant-Higgs, massless-neutrino current projection.
        coefficient = 0
        for node in restored.momentum_contributions[0].cut_side:
            f = leaves[node]
            f = f.strip_derivs() if f.derivs else f
            if f.label.rstrip("†") == "phidm":
                coefficient += -1 if f.is_conj else 1
        assert abs(coefficient) == expected
    assert not inspect_endpoints(raw[0]).agrees  # migration never mutates input
