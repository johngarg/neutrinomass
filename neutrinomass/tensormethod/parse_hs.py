#!/usr/bin/env python3

"""Tools for parsing the Hilbert Series into calls to ``invariants``."""

from copy import copy
from functools import reduce
from itertools import product

import sympy

import neutrinomass.tensormethod.hs as hs
import neutrinomass.tensormethod.sm as sm
import neutrinomass.tensormethod.core as tm
from neutrinomass.tensormethod.contract import invariants
from neutrinomass.tensormethod.hs import X

# plug in 3 fermion generations
H7_LNV_NF3 = hs.H7_LNV.xreplace({hs.Nf: 3})
H9_LNV_NF3 = hs.H9_LNV.xreplace({hs.Nf: 3})
H11_LNV_NF3 = hs.H11_LNV.xreplace({hs.Nf: 3})

FIELD_LOOKUP = {
    hs.L(X): sm.L,
    hs.Ld(X): sm.L.conj,
    hs.H(X): sm.H,
    hs.Hd(X): sm.H.conj,
    hs.Q(X): sm.Q,
    hs.Qd(X): sm.Q.conj,
    hs.eb(X): sm.eb,
    hs.ebd(X): sm.eb.conj,
    hs.ub(X): sm.ub,
    hs.ubd(X): sm.ub.conj,
    hs.db(X): sm.db,
    hs.dbd(X): sm.db.conj,
    hs.G(X): sm.G,
    hs.Gb(X): sm.Gb,
    hs.W(X): sm.W,
    hs.Wb(X): sm.Wb,
    hs.B(X): sm.B,
    hs.Bb(X): sm.Bb,
}


def distribute_derivatives(expr):
    """Returns a new Hilbert Series with the derivatives distributed across each
    term.

    For a single term, pass it in wrapped in a list.

    """
    new_terms = []
    f = lambda x: x.args if not isinstance(expr, list) else x
    for term in f(expr):
        # derivatives will never be outside of Mul
        if not isinstance(term, sympy.Mul):
            new_terms.append(term)
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


def is_number(expr):
    return isinstance(expr, sympy.Integer) or isinstance(expr, sympy.Rational)


def is_field(expr, term):
    if term:
        if isinstance(expr, sympy.Pow):
            expr = expr.args[0]

    return isinstance(expr, sympy.Function)


def is_deriv(expr):
    if isinstance(expr, sympy.Pow):
        expr = expr.args[0]

    return isinstance(expr, sympy.Derivative)


def is_term(expr):
    return isinstance(expr, sympy.Pow) or isinstance(expr, sympy.Mul)


def proc_number(expr):
    return [1]


def proc_field(expr):
    if isinstance(expr, sympy.Function):
        return [FIELD_LOOKUP[expr]]

    if isinstance(expr, sympy.Pow):
        base, power = expr.args
        return [FIELD_LOOKUP[base]] * power


def proc_deriv(expr):
    if isinstance(expr, sympy.Derivative):
        field, (_, n) = expr.args
        return [("D", n, FIELD_LOOKUP[field])]

    if isinstance(expr, sympy.Pow):
        base, power = expr.args
        return proc_deriv(base) * power


def is_sum(expr):
    return isinstance(expr, sympy.Add)


def is_symbolic_deriv(expr):
    # derivatives represented by tuples:
    # ("D", order, field)
    return isinstance(expr, tuple)


def no_numbers(expr):
    return [i for i in expr if not isinstance(i, int)]


def deriv_possibilities(field, order):
    if order < 1:
        return [field]

    if field.is_fermion:
        deltas = [(1, 1), (-1, 1), (1, -1)]
    else:
        deltas = [(1, 1), (-1, -1)]

    dynkin_options = []
    for delta_u, delta_d in deltas:
        u, d = field.lorentz_irrep
        sum_u = delta_u + u
        sum_d = delta_d + d
        if sum_u >= 0 and sum_d >= 0:
            u, d = delta_u + u, delta_d + d
            new_dynkin = str(u) + str(d)
            dynkin_options.append(new_dynkin)

    result = [deriv_possibilities(tm.D(field, d), order - 1) for d in dynkin_options]
    return sympy.flatten(result)


def proc_term(expr):
    flat_term = reduce(lambda x, y: x + y, expr)

    # expand derivative possibilities and find invariants, return as a list
    contains_deriv = False
    for item in flat_term:
        if is_symbolic_deriv(item):
            contains_deriv = True

    if not contains_deriv:
        # return [invariants(*no_numbers(flat_term))]
        return [no_numbers(flat_term)]

    # build new lists with derivative possibilities
    new_terms = [[]]
    for i, item in enumerate(flat_term):
        if not is_symbolic_deriv(item):
            for new_term in new_terms:
                new_term.append(item)
        if is_symbolic_deriv(item):
            _, order, field = item
            possible_fields = deriv_possibilities(field, order)
            new_terms = list(product(new_terms, possible_fields))
            # product leaves the list a bit dirty, need to clean:
            # ([ old list ], new_field) -> [*old_list, new_field]
            new_terms = [[*old_list, new_field] for old_list, new_field in new_terms]

    return [no_numbers(term) for term in new_terms]


def proc_sum(expr):
    return reduce(lambda x, y: x + y, expr)


def parse_hs(expr, term=False):
    if is_number(expr):
        return proc_number(expr)

    if is_field(expr, term=term):
        return proc_field(expr)

    if is_deriv(expr):
        return proc_deriv(expr)

    # recursive calls
    # term is a product of fields (not power)
    if is_term(expr):
        args = expr.args if not isinstance(expr, sympy.Pow) else [expr]
        return proc_term([parse_hs(item, term=True) for item in args])

    if is_sum(expr):
        return proc_sum([parse_hs(item) for item in expr.args])

    raise Exception(f"Missed a case for {expr} when parsing Hilbert Series.")


def parse(hs):
    """Parses Hilbert Series into a list of lists of fields."""
    return parse_hs(distribute_derivatives(hs))
