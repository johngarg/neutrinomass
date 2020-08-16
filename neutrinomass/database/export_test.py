#!/usr/bin/env python3

from neutrinomass.database.export import *
from neutrinomass.tensormethod import D, L, Q, db, ub, eb, H, eps, delta
from neutrinomass.completions.core import cons_completion_field
from neutrinomass.completions import operator_completions
from neutrinomass.completions import EFF_OPERATORS, DERIV_EFF_OPERATORS
from sympy import Rational

from networkx import Graph
from neutrinomass.completions import EffectiveOperator
from neutrinomass.database.database import LazyCompletion


ExoticField = cons_completion_field


def Partition(*args):
    return tuple(args)


def test_export_tensor():
    assert eval(export_tensor(L.conj("d0 i0"))).y == Rational(1, 2)

    # deriv example
    deriv_test_lhs = eval(export_tensor(D(L.conj, "10")("u0 i0")))
    deriv_test_rhs = D(L.conj, "10")("u0 i0")
    assert deriv_test_lhs == deriv_test_rhs

    # exotic example
    s1 = Field("S1", "00010", charges={"y": Rational(1, 3)})
    exotic_field = cons_completion_field(s1("-c0"))
    assert cons_completion_field(eval(export_tensor(exotic_field))) == exotic_field


def test_export_operator():
    op = (
        L("u0 i1")
        * Q("u1 c0 i2")
        * db("u2 -c1")
        * eps("-u0 -u1")
        * eps("-i1 -i2")
        * delta("c1 -c0")
    )
    assert eval(export_operator(op)) == op


def test_export_completion():
    for _, op in list(EFF_OPERATORS.items())[:5]:
        comp = list(operator_completions(op))[0]
        new_str, old = export_completion(comp, lazy=False), comp
        new = eval(new_str)

        assert new.exotics == old.exotics
        assert new.partition == old.partition
        assert new.graph.__dict__["_adj"] == old.graph.__dict__["_adj"]
        assert new.terms == old.terms
        assert new.operator.__dict__ == old.operator.__dict__


def test_lazy_completion():
    op = EFF_OPERATORS["3b"]
    comp = list(operator_completions(op))[0]
    lazy_comp = eval(export_completion(comp, lazy=True))

    assert isinstance(lazy_comp.head, dict)
    assert isinstance(lazy_comp.tail, str)

    op = DERIV_EFF_OPERATORS["D3"]
    comp = list(operator_completions(op))[0]
    lazy_comp = eval(export_completion(comp, lazy=True))
