#!/usr/bin/env python3

import sympy
from itertools import combinations, chain

from neutrinomass.database.closures import (
    neutrino_mass_estimate,
    numerical_np_scale_estimate,
    numerical_mv,
)


def get_leading_mv(eff_op):
    """Returns the symbolic expression for the leading order contribution to the
    neutrino mass implied by the operator.

    """
    estimates = neutrino_mass_estimate(eff_op)
    numerical = [(e, numerical_mv(e)) for e in estimates]
    return max(numerical, key=lambda x: x[1])[0]


def estimate_np_scale(eff_op):
    estimates = neutrino_mass_estimate(eff_op)
    numerical = [numerical_np_scale_estimate(e) for e in estimates]
    return numerical


def loop_data(expr):
    n_loops, n_loops_v2, max_loops = 0, [], 8
    for n in range(1, max_loops):
        if expr.coeff(sympy.Symbol("loop") ** n):
            n_loops += n

        if expr.coeff(sympy.Symbol("loopv2") ** n):
            n_loops_v2 += sorted([n - j for j in range(n + 1)])

    return n_loops, n_loops_v2


def table_data(eff_op):
    estimates = neutrino_mass_estimate(eff_op)
    numerical = [(e, numerical_np_scale_estimate(e)) for e in estimates]
    expr, np_scale = sorted(numerical, key=lambda x: round(x[1], 2))[-1]

    n_loops, n_loops_v2 = loop_data(expr)

    n_loops_str = (
        str(n_loops)
        if not n_loops_v2
        else ",".join(str(i + n_loops) for i in n_loops_v2)
    )
    return n_loops_str, f"\\mynum{{{np_scale}}}"


def subsets(lst):
    return list(chain(*[combinations(lst, i) for i in range(len(lst) + 1)]))
