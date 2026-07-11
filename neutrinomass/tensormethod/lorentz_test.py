#!/usr/bin/env python3

import pytest

from neutrinomass.completions.operators import DERIV_EFF_OPERATORS
from neutrinomass.tensormethod.contract import lorentz_singlets
from neutrinomass.tensormethod.core import Index, Operator
from neutrinomass.tensormethod.lorentz import (
    LorentzBasis,
    lorentz_contraction,
)


MULTIDIMENSIONAL_DERIVATIVE_OPERATORS = (
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
)


def renamed_lorentz_dummies(operator):
    replacements = []
    seen = set()
    for tensor in operator.tensors:
        for index in tensor.indices:
            if index.index_type not in Index.get_lorentz_index_types().values():
                continue
            raised = index if index.is_up else -index
            key = raised.index_type, str(raised)
            if key in seen:
                continue
            seen.add(key)
            short_type = (
                "u"
                if raised.index_type == Index.get_index_types()["u"]
                else "d"
            )
            fresh = Index.fresh(short_type)
            replacements.extend(((raised, fresh), (-raised, -fresh)))
    return Operator(
        *(tensor.fun_eval(*replacements) for tensor in operator.tensors)
    )


def test_d6a_lorentz_basis_has_rank_two_and_resolves_schouten_relation():
    operator = DERIV_EFF_OPERATORS["D6a"].operator
    basis = LorentzBasis.from_operator(operator)
    coordinates = {
        basis.project(singlet) for singlet in lorentz_singlets(operator)
    }

    assert basis.dimension == 2
    assert basis.labels == (
        "u:0-1,2-3|d:0-1",
        "u:0-2,1-3|d:0-1",
    )
    assert coordinates == {(1, 0), (0, 1), (1, -1)}


def test_ambient_port_basis_is_independent_of_field_statistics():
    basis = LorentzBasis.from_port_counts((("u", 4), ("d", 2)))

    assert basis.dimension == 2
    assert basis.labels == (
        "u:0-1,2-3|d:0-1",
        "u:0-2,1-3|d:0-1",
    )


@pytest.mark.parametrize("operator_name", MULTIDIMENSIONAL_DERIVATIVE_OPERATORS)
def test_every_multidimensional_operator_has_a_resolved_rank_two_lorentz_basis(
    operator_name,
):
    operator = DERIV_EFF_OPERATORS[operator_name].operator
    singlets = lorentz_singlets(operator)
    basis = LorentzBasis.from_operator(operator)

    assert len(singlets) == 3
    assert basis.dimension == 2
    assert all(basis.project(singlet) is not None for singlet in singlets)


def test_lorentz_projection_ignores_tensor_order_and_dummy_names():
    operator = DERIV_EFF_OPERATORS["D6a"].operator
    basis = LorentzBasis.from_operator(operator)
    singlet = lorentz_singlets(operator)[0]
    expected = basis.project(singlet)

    assert basis.project(Operator(*reversed(singlet.tensors))) == expected
    assert basis.project(renamed_lorentz_dummies(singlet)) == expected


def test_explicit_basis_component_has_an_orthogonal_negative_control():
    operator = DERIV_EFF_OPERATORS["D6a"].operator
    basis = LorentzBasis.from_operator(operator)
    by_label = {
        lorentz_contraction(singlet, normalise=True).label: singlet
        for singlet in lorentz_singlets(operator)
    }
    requested = basis.project(by_label[basis.labels[0]])
    orthogonal = basis.project(by_label[basis.labels[1]])

    assert requested[0] != 0
    assert orthogonal[0] == 0
