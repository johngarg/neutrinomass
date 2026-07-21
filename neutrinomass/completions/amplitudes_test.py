#!/usr/bin/env python3

from importlib import import_module
from pathlib import Path

import networkx as nx
import pytest
from sympy import Rational

from neutrinomass.completions.amplitudes import (
    audit_amplitude_symmetrisation,
    external_leg_provenance,
    induced_amplitude_operator,
)
from neutrinomass.completions.completions import (
    deriv_operator_completions,
    get_topology_data,
)
from neutrinomass.completions.operators import DERIV_EFF_OPERATORS
from neutrinomass.completions.core import (
    Completion,
    ComplexScalar,
    EffectiveOperator,
)
from neutrinomass.completions.topologies import Leaf
from neutrinomass.database import iter_completion_jsonl, write_completion_jsonl
from neutrinomass.database.rebuild import audit_amplitude_artifact
from neutrinomass.tensormethod.core import FERMI, IndexedField, Operator, eps


completions_module = import_module("neutrinomass.completions.completions")


JULIAN_MODEL = tuple(
    sorted(
        (
            ("F", 0, 0, 3, ("3b", 0), ("y", Rational("1/2"))),
            ("S", 0, 0, 2, ("3b", 0), ("y", 0)),
            ("S", 0, 0, 3, ("3b", 0), ("y", Rational("3/2"))),
        )
    )
)


VANISHING_MODEL = tuple(
    sorted(
        (
            ("F", 0, 0, 1, ("3b", 0), ("y", Rational("3/2"))),
            ("S", 0, 0, 2, ("3b", 0), ("y", 2)),
            ("S", 0, 0, 3, ("3b", 0), ("y", Rational("3/2"))),
        )
    )
)


@pytest.fixture(scope="module")
def d20_amplitude_witnesses():
    topology_data = [
        data
        for data in get_topology_data(5, 2)
        if Path(data["partition_file"]).stem in {"5s2f_9", "5s2f_11"}
    ]
    original = completions_module.get_topology_data
    completions_module.get_topology_data = lambda **kwargs: topology_data
    try:
        completions = deriv_operator_completions(DERIV_EFF_OPERATORS["D20"])
    finally:
        completions_module.get_topology_data = original

    julian = [
        completion
        for completion in completions
        if completion.topology == "5s2f_3"
        and tuple(sorted(completion.exotic_info().values())) == JULIAN_MODEL
    ]
    vanishing = [
        completion
        for completion in completions
        if completion.topology == "5s2f_24"
        and tuple(sorted(completion.exotic_info().values())) == VANISHING_MODEL
    ]
    assert julian
    assert vanishing
    return julian[0], vanishing[0]


def test_external_leg_provenance_is_complete_and_stable(d20_amplitude_witnesses):
    for completion in d20_amplitude_witnesses:
        provenance = external_leg_provenance(completion)
        assert len(provenance) == 7
        assert {leg.node for leg in provenance.values()} == {
            leaf.node
            for leaf in completions_module.partition_leaves(completion.partition)
        }


def test_external_leg_provenance_matches_repeated_index_free_fields_by_vertex():
    charges = {"y": 0, "3b": 0}
    eb_0 = IndexedField("eb†", "d0", charges=charges, comm=FERMI, is_conj=True)
    eb_1 = IndexedField("eb†", "d1", charges=charges, comm=FERMI, is_conj=True)
    ub_0 = IndexedField(
        "ub†", "d2 c0", charges=charges, comm=FERMI, is_conj=True
    )
    ub_1 = IndexedField(
        "ub†", "d3 c1", charges=charges, comm=FERMI, is_conj=True
    )
    phi = ComplexScalar("phi", "", charges=charges)
    graph = nx.Graph()
    graph.add_edge(1, 20, particle="eb†0")
    graph.add_edge(2, 20, particle="ub†0")
    graph.add_edge(3, 10, particle="eb†0")
    graph.add_edge(4, 10, particle="ub†0")
    graph.add_edge(10, 20, particle="phi")
    completion = Completion(
        operator=EffectiveOperator(
            "test", Operator(eb_0, ub_0, eb_1, ub_1)
        ),
        partition=(
            Leaf(eb_1, 1),
            Leaf(ub_1, 2),
            Leaf(eb_0, 3),
            Leaf(ub_0, 4),
        ),
        graph=graph,
        exotics={phi},
        terms=(
            Operator(phi.conj, eb_0, ub_0),
            Operator(phi, eb_1, ub_1),
        ),
        topology="test",
    )

    provenance = external_leg_provenance(completion)

    assert {leg.node for leg in provenance.values() if leg.term_index == 0} == {
        3,
        4,
    }
    assert {leg.node for leg in provenance.values() if leg.term_index == 1} == {
        1,
        2,
    }


def test_legacy_derivative_provenance_restores_generation_from_partition_leaf():
    charges = {"y": 0, "3b": 0}
    derivative = IndexedField(
        "Dx†",
        "u0 -c0 i0 g0",
        charges=charges,
        comm=FERMI,
        is_conj=True,
        derivs=1,
        stripped={
            "label": "x",
            "dynkin": "01011",
            "symmetry": [[1], [1], [1]],
            "charges": charges,
            "latex": "x",
        },
    )
    legacy_term_field = IndexedField(
        "x†",
        "d0 -c0 i0",
        charges=charges,
        comm=FERMI,
        is_conj=True,
    )
    companion = IndexedField(
        "y",
        "u1 c0 i1 g1",
        charges=charges,
        comm=FERMI,
    )
    legacy_companion = IndexedField(
        "y",
        "d1 c0 i1",
        charges=charges,
        comm=FERMI,
    )
    graph = nx.Graph()
    graph.add_edge(1, 10, particle="Dx†0")
    graph.add_edge(2, 10, particle="y0")
    completion = Completion(
        operator=EffectiveOperator(
            "test",
            Operator(
                derivative,
                companion,
                eps("-u0 -u1"),
                eps("-i0 -i1"),
            ),
        ),
        partition=(Leaf(derivative, 1), Leaf(companion, 2)),
        graph=graph,
        exotics=set(),
        terms=(
            Operator(
                legacy_term_field,
                legacy_companion,
                eps("-d0 -d1"),
                eps("-i0 -i1"),
            ),
        ),
        topology="test",
    )

    provenance = external_leg_provenance(completion)
    amplitude = induced_amplitude_operator(completion, {1: 1})
    restored = amplitude.indexed_fields[0]

    assert provenance[(0, 0)].node == 1
    assert restored.label == "Dx†"
    assert any(index.index_type == "Generation" for index in restored.indices)


def test_full_amplitude_symmetrisation_keeps_julian_and_rejects_counterexample(
    d20_amplitude_witnesses,
):
    julian, vanishing = d20_amplitude_witnesses

    julian_audit = audit_amplitude_symmetrisation(julian)
    vanishing_audit = audit_amplitude_symmetrisation(vanishing)

    assert all(term.safe_simplify() != 0 for term in vanishing.terms)
    assert julian_audit.status == "nonzero"
    assert julian_audit.witness
    assert vanishing_audit.status == "zero"


def test_amplitude_artifact_is_between_structural_and_physical_exact_classes(
    d20_amplitude_witnesses, tmp_path
):
    source = tmp_path / "structural.jsonl"
    destination = tmp_path / "physical.jsonl"
    julian, vanishing = d20_amplitude_witnesses
    write_completion_jsonl(source, [vanishing, julian])

    report = audit_amplitude_artifact(source, destination)
    survivors = list(iter_completion_jsonl(destination))

    assert report["input_records"] == 2
    assert report["rejected_records"] == 1
    assert report["surviving_records"] == 1
    assert [completion.topology for completion in survivors] == ["5s2f_3"]
