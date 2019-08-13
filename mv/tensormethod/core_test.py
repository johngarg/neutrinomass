#!/usr/bin/env python

from core import *

A = Field("A", dynkin="10011", charges={"y": 1})


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
        Field("AA", "11112", charges={"y": 0}),
        Field("AA", "11110", charges={"y": 0}),
        Field("AA", "11002", charges={"y": 0}),
        Field("AA", "11000", charges={"y": 0}),
    ]
