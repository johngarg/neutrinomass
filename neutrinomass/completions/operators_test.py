#!/usr/bin/env python3

from neutrinomass.completions.operators import EFF_OPERATORS, DERIV_EFF_OPERATORS


def test_pickle():
    op = DERIV_EFF_OPERATORS["D3"]
    assert sum(f.derivs for f in op.operator.fields)
