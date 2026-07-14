#!/usr/bin/env python3

from importlib import import_module
from pathlib import Path

import pytest
from sympy import Rational

from neutrinomass.completions.amplitudes import (
    audit_amplitude_symmetrisation,
    external_leg_provenance,
)
from neutrinomass.completions.completions import (
    deriv_operator_completions,
    get_topology_data,
)
from neutrinomass.completions.operators import DERIV_EFF_OPERATORS
from neutrinomass.database import iter_completion_jsonl, write_completion_jsonl
from neutrinomass.database.rebuild import audit_amplitude_artifact


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
