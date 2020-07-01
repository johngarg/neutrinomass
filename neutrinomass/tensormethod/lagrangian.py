#!/usr/bin/env python3

"""Provides the Lagrangian class and associated functions for symmetry
finding.

"""

from itertools import combinations_with_replacement
from neutrinomass.tensormethod.contract import invariants, unsimplified_invariants
from neutrinomass.completions.core import FieldType
from alive_progress import alive_bar
from multiset import Multiset
import numpy as np
from neutrinomass.tensormethod.core import (
    Operator,
    GENERATION,
    Index,
    decompose_product,
    Prod,
)

from sympy import Matrix
from sympy.tensor.tensor import tensorhead

# for now
from neutrinomass.tensormethod.sm import *

SM_MATTER = (H, Q, ub, db, L, eb)


def make_coupling(symbol: str, indices: str, sym=None):
    indices_list = indices.split()
    if sym is None:
        sym = [[1]] * len(indices_list)
    return tensorhead(symbol, [Index(i) for i in indices_list], sym)


def npoint_fieldstrings(n, fields=(L, eb, Q, db, ub, H), derivs=False, func=None):
    L.charges["l"] = 1
    eb.charges["l"] = -1
    Q.charges["l"] = 0
    ub.charges["l"] = 0
    db.charges["l"] = 0
    H.charges["l"] = 0

    conjs = tuple([f.conj for f in fields])
    if derivs:
        T = Field("D", "11000", charges={"y": 0, "3b": 0, "l": 0})
        fields += (T,)

    combos = list(combinations_with_replacement(fields + conjs, n))
    terms = []
    with alive_bar(len(combos)) as bar:
        while combos:
            combo = combos.pop(0)

            if func is not None:
                if not func(combo):
                    bar()
                    continue

            prods = decompose_product(*combo)
            singlets = [i for i in prods if i.is_singlet]
            if singlets:
                terms += [singlets[0]]
            bar()

    # for f in fields:
    #     del f.charges["l"]

    return terms


def prod_mass_dim(prod: Prod) -> int:
    """Returns the mass dimension of a Prod object"""
    out = []
    while True:
        irrep, left, right = prod
        out.append(right)

        if not isinstance(left, Prod):
            out.append(left)
            return sum(map(lambda x: x.mass_dim, out))

        prod = left


def generate_fieldstrings(max_dim, fields, derivs=True):
    """Returns a product of fields of maximum dimension `max_dim` built out of
    `fields`.

    """
    out = []
    for n_fields in range(2, max_dim + 1):
        for fieldstring in npoint_fieldstrings(n_fields, fields, derivs=derivs):
            if prod_mass_dim(fieldstring.walked()) <= max_dim:
                out.append(fieldstring)

    return out


def npoint_terms(n, fields, nf=3):
    conjs = [f.conj for f in fields]
    combos = combinations_with_replacement([*fields, *conjs], n)
    terms = []
    for combo in combos:
        invs = unsimplified_invariants(*combo, ignore=[])
        terms += invs

    return terms


class Lagrangian:
    def __init__(self, exotics, terms):
        self.exotics = exotics
        self.terms = terms

    def u1_symmetries(self):
        exotics = list(self.exotics)
        extra_0s = [0 for _ in range(len(exotics))]

        #            H  Q  ub db L eb
        yukawas = [[1, 1, 1, 0, 0, 0], [-1, 1, 0, 1, 0, 0], [-1, 0, 0, 0, 1, 1]]
        new_yukawas = [list(yuk) + extra_0s for yuk in yukawas]

        exotic_indices = {k: v for k, v in zip(exotics, range(0, len(exotics)))}

        matrix = []
        for term in self.terms:
            matrix += [term_to_row(term, self.exotics, exotic_indices)]

        matrix += new_yukawas
        X = Matrix(matrix)
        return X.nullspace()

    def is_conserving(self, symmetry):
        pass

    def generate_full(self):
        pass


def term_to_row(term, exotics, exotic_indices):
    # {H: 0, Q: 1, ub: 2, db: 3, L: 4, lb: 5}
    index_dict = dict(zip(SM_MATTER, range(6)))
    n_fields = len(exotics) + 6
    row = np.zeros(n_fields)
    for field, mult in Counter(term.fields).items():
        if not isinstance(field, FieldType):
            row[index_dict[field]] += mult if not field.is_conj else -mult
        else:
            row[6 + exotic_indices[field]] += mult if not field.is_conj else -mult

    return [int(i) for i in row]
