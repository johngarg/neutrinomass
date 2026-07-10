#!/usr/bin/env python3

from neutrinomass.completions.utils import fixedpoint, multiset_equality


def test_fixedpoint():
    assert fixedpoint(lambda value: min(value + 1, 3), 0) == 3


def test_multiset_equality():
    assert multiset_equality(["a", "a", "b"], ["b", "a", "a"])
    assert not multiset_equality(["a", "a", "b"], ["a", "b", "b"])
