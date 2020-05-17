#!/usr/bin/env python3

from neutrinomass.completions.core import *
from neutrinomass.completions.operators import EFF_OPERATORS, DERIV_EFF_OPERATORS
from neutrinomass.tensormethod.core import BOSE, FERMI, D
from neutrinomass.tensormethod.contract import invariants
from neutrinomass.tensormethod.sm import L, H


def test_complex_scalar():
    f = ComplexScalar("f", "i1 c1", charges={"y": 1})
    assert not f.is_conj
    assert f.y == 1
    assert f.comm == BOSE

    assert f.conj.is_conj
    assert f.conj.conj == f


def test_real_scalar():
    f = RealScalar("f", "i1")
    assert f.y == 0


def test_majorana_fermion():
    f = MajoranaFermion("f", "u0 c0 -c1")
    assert not f.is_conj
    assert f.conj.conj == f
    assert f.conj.is_conj
    assert invariants(f.field, f.field, ignore=[])


def test_vectorlike_dirac_fermion():
    f = VectorLikeDiracFermion("f", "u0 i1", charges={"y": 1})
    partner = f.dirac_partner()

    assert not f.is_conj
    assert f.conj.conj == f
    assert f.conj.is_conj
    assert invariants(f.field, f.dirac_partner().field, ignore=[])
    assert not invariants(f.field, f.field, ignore=[])
    return f, partner


def test_effective_operator():
    assert EFF_OPERATORS["1"].mass_dimension == 5
    assert EFF_OPERATORS["2"].mass_dimension == 7
    assert EFF_OPERATORS["11b"].mass_dimension == 9
    assert EFF_OPERATORS["47a"].mass_dimension == 11

    assert EFF_OPERATORS["11b"].topology_type == {"n_scalars": 0, "n_fermions": 6}
    assert EFF_OPERATORS["1"].topology_type == {"n_scalars": 2, "n_fermions": 2}

    assert EFF_OPERATORS["1"].indexed_fields
    assert EFF_OPERATORS["11a"].fields


def test_topology_type():
    """(DL)^{i} L^{j} (DH)^{k} (DDH)^{l}"""
    eff_ops = [
        EffectiveOperator("", op)
        for op in invariants(D(L, "01"), L, D(H, "11"), D(D(H, "11"), "22"))
    ]
    assert eff_ops[0].topology_type == {"n_scalars": 2, "n_fermions": 2}
