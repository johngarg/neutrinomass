#!/usr/bin/env python3

"""Provides the Lagrangian class and associated functions for symmetry
finding.

"""

from itertools import combinations_with_replacement
from neutrinomass.tensormethod.contract import invariants
from alive_progress import alive_bar
from neutrinomass.tensormethod.core import (
    Operator,
    GENERATION,
    Index,
    decompose_product,
)

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


def npoint_fieldstrings(n, fields=[L, eb, Q, db, ub, H], derivs=True):
    if derivs:
        T = Field("D", "11000")
        fields += [T]

    conjs = [f.conj for f in fields]
    combos = list(combinations_with_replacement(fields + conjs, n))
    terms = []
    with alive_bar(len(combos)) as bar:
        for combo in combos:
            prods = decompose_product(*combo)
            singlets = [i for i in prods if i.is_singlet]
            if singlets:
                terms.append(singlets[0])
            bar()

    return terms


def npoint_terms(n, fields, nf=3):
    conjs = [f.conj for f in fields]
    combos = combinations_with_replacement(fields + conjs, n)
    terms = []
    for combo in combos:
        invs = invariants(*combo, ignore=["c", "u", "d"])
        terms += invs

    return terms


class Lagrangian:
    def __init__(self, *terms):
        self.terms = terms
