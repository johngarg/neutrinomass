#!/usr/bin/env python3

from neutrinomass.completions.completion import *
from neutrinomass.completions.operators import EFF_OPERATORS


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
    o2_comps = filter_completions(o2, sieve)
    o3a_comps = filter_completions(o3a, sieve)
    o3b_comps = filter_completions(o3b, sieve)
    o4a_comps = filter_completions(o4a, sieve)
    o4b_comps = filter_completions(o4b, sieve)
    o8_comps = filter_completions(o8, sieve)

    assert len(o2_comps) == 2
    assert len(o3a_comps) == 6
    assert len(o3b_comps) == 5
    assert len(o4a_comps) == 0
    assert len(o4b_comps) == 3
    assert len(o8_comps) == 4
