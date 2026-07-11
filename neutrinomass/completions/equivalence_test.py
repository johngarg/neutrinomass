#!/usr/bin/env python3

import networkx as nx
from sympy import Rational

from neutrinomass.completions.completions import (
    are_equivalent_completions,
    check_remapping_on_terms,
    compare_terms,
)
from neutrinomass.completions.core import (
    Completion,
    ComplexScalar,
    EffectiveOperator,
    RealScalar,
    VectorLikeDiracFermion,
)
from neutrinomass.completions.equivalence import equivalent_lagrangians
from neutrinomass.completions.operators import EFF_OPERATORS
from neutrinomass.tensormethod import H, L, eps
from neutrinomass.tensormethod.core import IndexedField


def make_completion(exotics, terms):
    return Completion(
        operator=EffectiveOperator("test", EFF_OPERATORS["1"].operator),
        partition=(),
        graph=nx.Graph(),
        exotics=set(exotics),
        terms=list(terms),
        topology="test_1",
    )


def doublet(label, index):
    return ComplexScalar(
        label,
        index,
        charges={"y": Rational(1, 2), "3b": 0},
    )


def higgs_interaction(field, higgs, conjugate_higgs=False):
    higgs_field = H.conj(higgs) if conjugate_higgs else H(higgs)
    return field.conj * higgs_field * eps(f"-{field.indices[0]} -{higgs}")


def test_equivalence_enumerates_identical_representation_bijections():
    a, b = doublet("a", "i0"), doublet("b", "i1")
    x, y = doublet("x", "i2"), doublet("y", "i3")
    m, n = doublet("m", "i4"), doublet("n", "i5")
    first = make_completion(
        [a, b],
        [higgs_interaction(a, "i10"), higgs_interaction(b, "i11", True)],
    )
    second = make_completion(
        [x, y],
        [higgs_interaction(y, "i12"), higgs_interaction(x, "i13", True)],
    )
    third = make_completion(
        [m, n],
        [higgs_interaction(n, "i14"), higgs_interaction(m, "i15", True)],
    )

    assert compare_terms(first, second) == {"a": "y", "b": "x"}
    assert are_equivalent_completions(first, first)
    assert are_equivalent_completions(first, second)
    assert are_equivalent_completions(second, first)
    assert are_equivalent_completions(second, third)
    assert are_equivalent_completions(first, third)


def test_equivalence_preserves_distinct_species_multiplicity():
    a, b = doublet("a", "i0"), doublet("b", "i1")
    x = doublet("x", "i2")
    two_species = make_completion(
        [a, b],
        [higgs_interaction(a, "i10"), higgs_interaction(b, "i11")],
    )
    one_species = make_completion(
        [x],
        [higgs_interaction(x, "i12"), higgs_interaction(x, "i13")],
    )

    assert not are_equivalent_completions(two_species, one_species)


def test_term_comparison_preserves_inequivalent_su2_contractions():
    a = IndexedField("a", "i0", charges={"y": 0, "3b": 0})
    b = IndexedField("b", "i1", charges={"y": 0, "3b": 0})
    c = IndexedField("c", "i2", charges={"y": 0, "3b": 0})
    d = IndexedField("d", "i3", charges={"y": 0, "3b": 0})
    first = a * b * c * d * eps("-i0 -i1") * eps("-i2 -i3")
    second = a * b * c * d * eps("-i0 -i2") * eps("-i1 -i3")

    assert check_remapping_on_terms(
        [first], [second], {"a": "a", "b": "b", "c": "c", "d": "d"}
    ) == {}


def test_term_comparison_preserves_interaction_multiplicity():
    a, b = doublet("a", "i0"), doublet("b", "i1")
    first = higgs_interaction(a, "i10")
    second = higgs_interaction(b, "i11")

    assert equivalent_lagrangians([first, first], [first, first])
    assert not equivalent_lagrangians([first, first], [first, second])


def test_remapping_matches_whole_labels_not_prefixes():
    a = doublet("a", "i0")
    a0 = doublet("a0", "i1")
    x = doublet("x", "i2")
    y = doublet("y", "i3")

    assert check_remapping_on_terms(
        [higgs_interaction(a, "i10"), higgs_interaction(a0, "i11")],
        [higgs_interaction(x, "i12"), higgs_interaction(y, "i13")],
        {"a": "x", "a0": "y"},
    ) == {"a": "x", "a0": "y"}


def test_equivalence_allows_conjugate_field_naming_conventions():
    positive = doublet("a", "i0")
    negative = ComplexScalar(
        "x",
        "i1",
        charges={"y": Rational(-1, 2), "3b": 0},
    )
    first = make_completion(
        [positive],
        [higgs_interaction(positive, "i10")],
    )
    second = make_completion(
        [negative],
        [negative * H("i11") * eps("-i1 -i11")],
    )

    assert compare_terms(first, second) == {"a": "x"}
    assert compare_terms(second, first) == {"x": "a"}


def test_equivalence_allows_dirac_partner_naming_conventions():
    source = VectorLikeDiracFermion(
        "a", "u0", charges={"y": 0, "3b": 0}
    )
    target = VectorLikeDiracFermion(
        "x", "u1", charges={"y": 0, "3b": 0}
    ).dirac_partner()
    first = make_completion([source], [source * H("i0")])
    second = make_completion([target], [target * H("i1")])

    assert compare_terms(first, second) == {"a": "x"}
    assert compare_terms(second, first) == {"x": "a"}


def test_equivalence_canonicalises_dirac_conjugate_suffix_order():
    field = VectorLikeDiracFermion(
        "a", "u0", charges={"y": 0, "3b": 0}
    )
    partner_then_conjugate = field.dirac_partner().conj
    conjugate_then_partner = field.conj.dirac_partner()

    assert partner_then_conjugate.label == "a~†"
    assert conjugate_then_partner.label == "a†~"
    assert equivalent_lagrangians(
        [partner_then_conjugate * H("i0")],
        [conjugate_then_partner * H("i1")],
    )


def test_equivalence_distinguishes_real_and_complex_particles():
    real = RealScalar("r", "i0", charges={"y": 0, "3b": 0})
    complex_ = ComplexScalar("c", "i1", charges={"y": 0, "3b": 0})
    real_completion = make_completion([real], [real * H("i2")])
    complex_completion = make_completion([complex_], [complex_ * H("i3")])

    assert not are_equivalent_completions(real_completion, complex_completion)


def test_graph_cache_preserves_custom_field_metadata():
    neutral = IndexedField("a", "i0", charges={"y": 0, "3b": 0})
    charged = IndexedField("a", "i0", charges={"y": 1, "3b": 0})
    neutral_term = neutral * H("i1") * eps("-i0 -i1")
    charged_term = charged * H("i1") * eps("-i0 -i1")

    assert not equivalent_lagrangians([neutral_term], [charged_term])


def test_equivalence_ignores_free_generation_labels():
    with_generation = (
        L("u0 i0 g0") * H("i1") * eps("-i0 -i1")
    )
    without_generation = L("u1 i2") * H("i3") * eps("-i2 -i3")

    assert equivalent_lagrangians([with_generation], [without_generation])


def test_interactions_are_equivalent_to_their_hermitian_conjugates():
    left = IndexedField("a", "u0", charges={"y": 1, "3b": 0})
    right = IndexedField("b", "u1", charges={"y": -1, "3b": 0})
    interaction = left * right * eps("-u0 -u1")
    conjugate = left.conj * right.conj * eps("-d0 -d1")

    assert equivalent_lagrangians([interaction], [conjugate])
