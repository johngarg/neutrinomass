#!/usr/bin/env python3

from mv.completions.completion import *
from mv.completions.operators import EFF_OPERATORS


def lnv_completions(op):
    return operator_completions(EFF_OPERATORS[op])


# completions of dimension 5 and 7 operators
sieve = collect_completions(lnv_completions("1"))
o2 = collect_completions(lnv_completions("2"))
o3a = collect_completions(lnv_completions("3a"))
o3b = collect_completions(lnv_completions("3b"))
o4a = collect_completions(lnv_completions("4a"))
o4b = collect_completions(lnv_completions("4b"))
o8 = collect_completions(lnv_completions("8"))


def test_completions():
    # test on results of 1410.0689
    assert len(filter_completions(o2, sieve)) == 2
    assert len(filter_completions(o3a, sieve)) == 6
    assert len(filter_completions(o3b, sieve)) == 5
    assert len(filter_completions(o4a, sieve)) == 0
    assert len(filter_completions(o4b, sieve)) == 3
    assert len(filter_completions(o8, sieve)) == 4
