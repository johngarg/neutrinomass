#!/usr/bin/env python3

from neutrinomass.tensormethod.lagrangian import *

from neutrinomass.tensormethod import eps
from neutrinomass.completions import EFF_OPERATORS, operator_completions
from neutrinomass.completions.core import cons_completion_field
from sympy import Rational

TEST_COMP = list(operator_completions(EFF_OPERATORS["1"]))[0]
LAG = Lagrangian(TEST_COMP.exotics, TEST_COMP.terms)


def test_u1_symmetries():
    exotics = [f.field for f in TEST_COMP.exotics]
    exotic_indices = {k: v for k, v in zip(exotics, range(0, len(exotics)))}
    exotic_indices = {
        **exotic_indices,
        **{k: v for k, v in zip([f.conj for f in exotics], range(0, len(exotics)))},
    }
    assert term_to_row(TEST_COMP.terms[0], TEST_COMP.exotics, exotic_indices)
    assert LAG.num_u1_symmetries() == 2


def test_generate_uv_terms():
    S1 = Field("S1", "00010", charges={"y": Rational(1, 3)})("-c0")
    s1_complex_scalar = cons_completion_field(S1)
    lag = Lagrangian(
        {s1_complex_scalar},
        [
            s1_complex_scalar
            * Q("u0 c0 i0")
            * L("u1 i1")
            * eps("-u0 -u1")
            * eps("-i0 -i1"),
            s1_complex_scalar * ub.conj("d0 c0") * eb.conj("d1") * eps("-d0 -d1"),
        ],
    )
    full_lag = lag.generate_full()

    assert lag.num_u1_symmetries() == 3
    assert full_lag.num_u1_symmetries() == 2
