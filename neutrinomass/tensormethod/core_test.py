#!/usr/bin/env python

from neutrinomass.tensormethod.core import *

A = Field("A", dynkin="10011", charges={"y": 1})
B = IndexedField("B", indices="i0 i1")
C = IndexedField("C", indices="c0 c1")
dd = IndexedField("D", indices="-c0 c2", charges={"y": 2, "3b": 1})


def test_get_dynkin():
    assert get_dynkin("u0 c0 -c1 i0") == "10111"
    assert get_dynkin(B.indices) == "00002"


def test_conj():
    assert A.conj.conj == A
    assert not A.conj == A
    assert A.conj.y == -1
    assert A.conj.conj.y != A.conj.y


def test_indices():
    assert type(A("u0 -c1 i0")) == IndexedField
    cnj_A = IndexedField("A", "d0 c1 i0", charges={"y": -1}, is_conj=True)
    assert A("u0 -c1 i0").conj.quantum_numbers == cnj_A.quantum_numbers
    assert A("u0 -c1 i0").conj.indices == cnj_A.indices


def test_mul():
    assert A * A.conj == [
        Field("AA†", "11112", charges={"y": 0}),
        Field("AA†", "11110", charges={"y": 0}),
        Field("AA†", "11002", charges={"y": 0}),
        Field("AA†", "11000", charges={"y": 0}),
    ]


def test_op():
    prod = C * dd
    assert len(prod.free_indices) == 2
    assert prod.dynkin == "00200"


def test_fresh():
    assert isinstance(A.fresh_indices(), IndexedField)


def test_decompose_product():
    prods = decompose_product(A, B.field, A.conj)

    first_two = [
        {
            "label": "ABA†",
            "dynkin": "11114",
            "charges": {"y": 0},
            "is_conj": False,
            "comm": 0,
            "symmetry": [[1], [1], [1], [1], [1, 1, 1, 1]],
        },
        {
            "label": "ABA†",
            "dynkin": "11112",
            "charges": {"y": 0},
            "is_conj": False,
            "comm": 0,
            "symmetry": [[1], [1], [1], [1], [1, 1]],
        },
    ]

    last_two = [
        {
            "label": "ABA†",
            "dynkin": "11002",
            "charges": {"y": 0},
            "is_conj": False,
            "comm": 0,
            "symmetry": [[1], [1], [1, 1]],
        },
        {
            "label": "ABA†",
            "dynkin": "11000",
            "charges": {"y": 0},
            "is_conj": False,
            "comm": 0,
            "symmetry": [[1], [1]],
        },
    ]
    assert prods[:2] == [Field(**x) for x in first_two]
    assert prods[-2:] == [Field(**x) for x in last_two]


def test_indexed_field_conj():
    assert B.conj.conj == B
    assert not B.is_conj
    assert B.conj.is_conj

    assert C.conj.conj == C
    assert not C.is_conj
    assert C.conj.is_conj

    assert dd.conj.conj == dd
    assert not dd.is_conj
    assert dd.conj.is_conj
    assert dd.conj.y == -2


def test_latex_preserves_direct_colour_contractions_and_variance():
    upper = IndexedField("A", "c0", latex="A")
    lower = IndexedField("B", "-c0", latex="B")

    assert (upper * lower).latex() == r"B_{a} A^{a}"


def test_latex_preserves_delta_and_colour_epsilon_variance():
    upper = IndexedField("A", "c0", latex="A")
    lower = IndexedField("B", "-c1", latex="B")
    with_delta = upper * lower * delta("c1 -c0")

    assert with_delta.latex() == r"B_{a} A^{b}  \cdot  \delta^{a}_{b}"

    three_upper = IndexedField("A", "c0 c1 c2", latex="A")
    three_lower = IndexedField("A", "-c0 -c1 -c2", latex="A")
    assert (three_upper * eps("-c0 -c1 -c2")).latex() == (
        r"A^{a b c}  \cdot  \epsilon_{a b c}"
    )
    assert (three_lower * eps("c0 c1 c2")).latex() == (
        r"A_{a b c}  \cdot  \epsilon^{a b c}"
    )


def test_latex_preserves_antisymmetric_index_order():
    first = IndexedField("A", "u0", latex="A")
    second = IndexedField("B", "u1", latex="B")
    ordered = first * second * eps("-u0 -u1")

    first = IndexedField("A", "u0", latex="A")
    second = IndexedField("B", "u1", latex="B")
    reversed_ = first * second * eps("-u1 -u0")

    assert ordered.latex() == (
        r"B^{\alpha} A^{\beta}  \cdot  \epsilon_{\beta \alpha}"
    )
    assert reversed_.latex() == (
        r"B^{\alpha} A^{\beta}  \cdot  \epsilon_{\alpha \beta}"
    )
    assert ordered.latex() != reversed_.latex()


def test_latex_ignore_skips_invariants_of_ignored_type():
    upper = IndexedField("A", "c0", latex="A")
    lower = IndexedField("B", "-c1", latex="B")
    assert (upper * lower * delta("c1 -c0")).latex(ignore="c") == "B A"

    first = IndexedField("A", "u0", latex="A")
    second = IndexedField("B", "u1", latex="B")
    assert (first * second * eps("-u0 -u1")).latex(ignore="u") == "B A"


def test_latex_displays_generation_indices_by_default():
    flavoured = IndexedField("A", "u0 g0", latex="A", nf=3)
    singlet = IndexedField("B", "", latex="B")

    assert (flavoured * singlet).latex() == r"B A^{\alpha p}"
    assert (flavoured * singlet).latex(ignore="g") == r"B A^{\alpha}"


def test_strip_derivs():
    from neutrinomass.tensormethod.sm import Q, H
    from neutrinomass.tensormethod.core import D

    stripped_q = D(D(Q, "01"), "10")(
        "u0 c0 i0 g0"
    ).strip_derivs_with_indices()
    expected_q = Q("u0 c0 i0 g0")
    assert stripped_q.label == expected_q.label
    assert stripped_q.indices == expected_q.indices
    assert stripped_q.charges == expected_q.charges
    assert stripped_q.derivs == 0

    stripped_h = D(D(H, "11"), "00")("i0").strip_derivs_with_indices()
    expected_h = H("i0")
    assert stripped_h.label == expected_h.label
    assert stripped_h.indices == expected_h.indices
    assert stripped_h.charges == expected_h.charges
    assert stripped_h.derivs == 0

    assert D(Q, "01")("d0 c0 i0").strip_derivs_with_indices().derivs == 0
