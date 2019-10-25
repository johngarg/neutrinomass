#!/usr/bin/env python

from core import *

A = Field("A", dynkin="10011", charges={"y": 1})
B = IndexedField("B", indices="i0 i1")
C = IndexedField("C", indices="c0 c1")
D = IndexedField("D", indices="-c0 c2")


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
    assert A("u0 -c1 i0").conj == cnj_A


def test_mul():
    assert A * A.conj == [
        Field("AA†", "11112", charges={"y": 0}),
        Field("AA†", "11110", charges={"y": 0}),
        Field("AA†", "11002", charges={"y": 0}),
        Field("AA†", "11000", charges={"y": 0}),
    ]


def test_op():
    prod = C * D
    assert len(prod.free_indices) == 2
    assert prod.dynkin == "00200"


def test_fresh():
    assert isinstance(A.fresh_indexed_field(), IndexedField)


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


# calls
test_get_dynkin()
test_conj()
test_indices()
test_mul()
test_op()
test_fresh()
test_decompose_product()
