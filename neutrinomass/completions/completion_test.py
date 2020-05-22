#!/usr/bin/env python3

from neutrinomass.completions.completion import *
from neutrinomass.completions.operators import EFF_OPERATORS


def lnv_completions(op):
    return operator_completions(EFF_OPERATORS[op])


SIEVE = collect_completions(lnv_completions("1"))


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
    assert len(o4a_comps) == 0
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


def test_deriv_lnv_completions():
    # Derivative operator examples:
    from neutrinomass.tensormethod.sm import L, db, H
    from neutrinomass.tensormethod.core import D

    od2 = EffectiveOperator(
        "OD2",
        L("u0 i0")
        * L("u1 i1")
        * H("i2")
        * D(H, "11")("u2 d0 i3")
        * db.conj("d3 c0")
        * db("u4 -c1")
        * eps("-i0 -i1")
        * eps("-i2 -i3"),
    )

    comps = collect_completions(operator_completions(od2))
    assert len(comps) == 22


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
