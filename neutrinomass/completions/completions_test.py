#!/usr/bin/env python3

from neutrinomass.completions.completions import *
from neutrinomass.tensormethod.sm import L, Q, H
from neutrinomass.tensormethod.core import D

from neutrinomass.completions.operators import EFF_OPERATORS, DERIV_EFF_OPERATORS

from sympy import Rational

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


def test_construct_completion():
    data = partitions(EFF_OPERATORS["2"])[0]
    assert construct_completion(data["partition"], data["epsilons"], data["graph"])


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
    for comps, models in [
        (o7_comps, o7_models_from_paper),
        (o9_comps, o9_models_from_paper),
    ]:
        for k, v in comps.items():
            if frozenset(k) not in models:
                assert False

        for model in models:
            if model not in set(map(frozenset, comps.keys())):
                assert False

    assert True
