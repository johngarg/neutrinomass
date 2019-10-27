#!/usr/bin/env python3

"""Tools for parsing the Hilbert Series into calls to ``tensormethod.invariants``."""

from hs import X
import hs
import sm
import sympy

from copy import copy

# plug in nf = 3
H7_LNV_NF3 = hs.H7_LNV.xreplace({hs.Nf: 3})
H9_LNV_NF3 = hs.H9_LNV.xreplace({hs.Nf: 3})


def distribute_derivatives(expr):
    new_terms = []
    for i, term in enumerate(expr.args):
        # derivatives will never be outside of Mul
        if not isinstance(term, sympy.Mul):
            continue
        # iterate through items in a term to extract derivative order if present
        for item in term.args:
            if not str(item).startswith("D"):
                continue

            if isinstance(item, sympy.Pow):
                base, power = item.args
                if base == hs.D:
                    new_term = term / (hs.D ** power)
                    for _ in range(power):
                        new_term = new_term.diff(X)
            else:
                new_term = term / hs.D
                new_term = new_term.diff(X)

            new_terms.append(new_term)
            break

        else:
            new_terms.append(term)

    return sum(new_terms)


def parse_hs(expr):
    pass


def generate_field_collections(hs):
    seen = set()
    out = []
    for field_coll in parse_hs(hs):
        # iterate over two outputs of apply_eom
        hashable = sorted(tuple(field_coll))
        if hashable not in seen:
            seen.add(hashable)
            out.append(field_coll)

    return out
