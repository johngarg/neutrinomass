from copy import copy
import sqlite3

from neutrinomass.completions.core import EffectiveOperator
from neutrinomass.completions.equivalence_test import (
    doublet, higgs_interaction, make_completion,
)
from neutrinomass.database.global_equivalence import (
    equivalent_global_lagrangians, GlobalLagrangianIndex,
)


def test_global_equivalence_ignores_operator_and_topology_but_preserves_terms():
    a, x = doublet("a", "i0"), doublet("x", "i1")
    left = make_completion([a], [higgs_interaction(a, "i10")])
    right = make_completion([x], [higgs_interaction(x, "i11")])
    right.operator = EffectiveOperator("different", right.operator.operator)
    right.topology = "another_topology"
    assert equivalent_global_lagrangians(left, right)
    repeated = copy(right)
    repeated.terms = right.terms * 2
    assert not equivalent_global_lagrangians(left, repeated)


def test_global_equivalence_preserves_copies():
    a, b = doublet("a", "i0"), doublet("b", "i1")
    x = doublet("x", "i2")
    two = make_completion([a,b], [higgs_interaction(a,"i10"), higgs_interaction(b,"i11")])
    one = make_completion([x], [higgs_interaction(x,"i12"), higgs_interaction(x,"i13")])
    assert not equivalent_global_lagrangians(one, two)


def test_global_comparison_preserves_contractions_and_conjugate_naming():
    from neutrinomass.tensormethod.core import IndexedField,eps
    from neutrinomass.tensormethod import H
    from neutrinomass.completions.core import ComplexScalar
    from sympy import Rational
    a,b,c,d=[IndexedField(label,f"i{i}",charges={"y":0,"3b":0})
             for i,label in enumerate("abcd")]
    first=a*b*c*d*eps("-i0 -i1")*eps("-i2 -i3")
    second=a*b*c*d*eps("-i0 -i2")*eps("-i1 -i3")
    assert not equivalent_global_lagrangians(make_completion([],[first]),make_completion([],[second]))
    positive=doublet("positive","i4")
    negative=ComplexScalar("negative","i5",charges={"y":Rational(-1,2),"3b":0})
    assert equivalent_global_lagrangians(
        make_completion([positive],[higgs_interaction(positive,"i6")]),
        make_completion([negative],[negative*H("i7")*eps("-i5 -i7")]))


def test_global_index_merges_relabelled_class_and_unions_filtered_status():
    from neutrinomass.completions.completions import operator_completions
    from neutrinomass.completions.operators import EFF_OPERATORS
    completion = next(operator_completions(EFF_OPERATORS["1"]))
    with sqlite3.connect(":memory:") as connection:
        index = GlobalLagrangianIndex(connection)
        first, added = index.add(completion)
        assert added
        second, added = index.add(completion, filtered=True)
        assert first == second and not added
        assert index.counts() == (1,1)
