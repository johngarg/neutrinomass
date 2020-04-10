#!/usr/bin/env python3

from neutrinomass.tensormethod.core import Index, delta, eps

from neutrinomass.tensormethod.contract import colour_singlets
from neutrinomass.tensormethod.contract import construct_operators
from neutrinomass.tensormethod.contract import unsimplified_invariants
from neutrinomass.tensormethod.contract import extract_relabellings
from neutrinomass.tensormethod.contract import compare_singlets
from neutrinomass.tensormethod.contract import is_identity_mapping
from neutrinomass.tensormethod.contract import clean_operators
from neutrinomass.tensormethod.contract import remove_relabellings
from neutrinomass.tensormethod.contract import invariants
from neutrinomass.tensormethod.sm import *


def test_colour_singlets():
    singlets = colour_singlets([Q("u0 c0 i0") * db("u1 -c1")])
    assert singlets[0] == Q("u0 c0 i0") * db("u1 -c1") * delta("c1 -c0")

    singlets = colour_singlets([G("u0 u1 c0 -c1") * G("u2 u3 c2 -c3")])
    assert len(singlets) == 2
    # tensormethod doesn't know that G should be traceless on c0 and c1
    ans1 = G("u0 u1 c0 -c1") * G("u2 u3 c2 -c3") * delta("c1 -c0") * delta("c3 -c2")
    ans2 = G("u0 u1 c0 -c1") * G("u2 u3 c2 -c3") * delta("c1 -c2") * delta("c3 -c0")
    assert singlets == [ans1, ans2]

    singlets = colour_singlets(
        [Q("u0 c0 i0") * db("u1 -c1"), Q("u0 c0 i0") * ub("u1 -c1")]
    )
    assert len(singlets) == 2


def test_construct_operators():
    prods = [i.walked() for i in L * L]
    ans = [
        [L("u0 i0") * L("u1 i1")],
        [L("u0 i0") * L("u1 i1") * eps("-i0 -i1")],
        [L("u0 i0") * L("u1 i1") * eps("-u0 -u1")],
        [L("u0 i0") * L("u1 i1") * eps("-u0 -u1") * eps("-i0 -i1")],
    ]

    # need to define a better operator equality to test this better
    assert len([construct_operators(i) for i in prods]) == len(ans)


def test_unsimplified_invariants():
    pass


def test_extract_relabellings():
    a = [eps("-i0 -i2"), eps("-i3 -i1")]
    b = [eps("-i3 -i2"), eps("-i2 -i1")]
    # relabellings action on a
    #    0 -> 2, 2 -> 3, 3 -> 2
    #    1 -> 2, 0 -> 1
    relabellings1 = [
        (Index("-i0"), Index("-i2")),
        (Index("-i2"), Index("-i3")),
        (Index("-i3"), Index("-i2")),
    ]
    relabellings2 = [(Index("-i1"), Index("-i2")), (Index("-i0"), Index("-i1"))]
    assert extract_relabellings(a, b) == [relabellings1, relabellings2]


def test_compare_singlets():
    pass


def test_is_identity_mapping():
    pass


def test_clean_operators():
    pass


def test_remove_relabellings():
    pass


def test_invariants():
    o1 = invariants(L, L, H, H)
    o2 = invariants(L, L, L, eb, H)
    o3 = invariants(L, L, Q, db, H)
    o4 = invariants(L, L, Q.conj, ub.conj, H)
    o5 = invariants(L, L, Q, db, H, H, H.conj)
    # o29 = invariants(L, L, Q, Q.conj, ub, ub.conj, H, H)
    assert len(o1) == 1
    assert len(o2) == 1
    assert len(o3) == 2
    assert len(o4) == 2
    assert len(o5) == 4
    # assert len(o29) == 4
