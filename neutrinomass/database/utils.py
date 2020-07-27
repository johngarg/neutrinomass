#!/usr/bin/env python3

import sympy
from neutrinomass.database.closures import (
    neutrino_mass_estimate,
    numerical_np_scale_estimate,
)


def get_leading_mv(eff_op):
    """Returns the symbolic expression for the leading order contribution to the
    neutrino mass implied by the operator.

    """
    estimates = neutrino_mass_estimate(eff_op)
    numerical = [(e, numerical_np_scale_estimate(e)) for e in estimates]
    return sorted(numerical, key=lambda x: x[1])[-1][0]


def estimate_np_scale(eff_op):
    estimates = neutrino_mass_estimate(eff_op)
    numerical = [numerical_np_scale_estimate(e) for e in estimates]
    return numerical


def table_data(eff_op):
    estimates = neutrino_mass_estimate(eff_op)
    numerical = [(e, numerical_np_scale_estimate(e)) for e in estimates]
    expr, np_scale = sorted(numerical, key=lambda x: x[1])[-1]

    expr.replace(sympy.Symbol("loopv2"), sympy.Symbol("loop"))
    n_loops, max_loops = 0, 8
    for n in range(1, max_loops):
        if expr.coeff(sympy.Symbol("loop") ** n):
            n_loops = n
            break

    return n_loops, f"\\num{{{np_scale}}}"
