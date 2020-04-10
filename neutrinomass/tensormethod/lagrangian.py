#!/usr/bin/env python3

"""Provides the Lagrangian class and associated functions for symmetry
finding.

"""

from itertools import combinations_with_replacement
from neutrinomass.tensormethod.contract import invariants
from neutrinomass.tensormethod.core import Operator, GENERATION, Index

from sympy.tensor.tensor import tensorhead

# for now
from neutrinomass.tensormethod.sm import *


def make_coupling(symbol: str, indices: str, sym=None):
    indices_list = indices.split()
    if sym is None:
        sym = [[1]] * len(indices_list)
    return tensorhead(symbol, [Index(i) for i in indices_list], sym)


def is_good_term(term: Operator):
    pass


def npoint_terms(n, fields, nf=3):
    conjs = [f.conj for f in fields]
    combos = combinations_with_replacement(fields + conjs, n)
    terms = []
    for combo in combos:
        invs = invariants(*combo, ignore=[])
        terms += invs

    return terms


class Lagrangian:
    def __init__(self, *terms):
        self.terms = terms
