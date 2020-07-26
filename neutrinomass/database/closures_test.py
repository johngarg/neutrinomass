#!/usr/bin/env python3

from neutrinomass.database.closures import *
from neutrinomass.completions import EFF_OPERATORS, DERIV_EFF_OPERATORS
from neutrinomass.database.utils import get_leading_mv, estimate_np_scale

loop = sympy.Symbol("loop")
loopv2 = sympy.Symbol("loopv2")
v = sympy.Symbol("v")
Λ = sympy.Symbol("Λ")
yd = sympy.Symbol("yd")
ye = sympy.Symbol("ye")
yu = sympy.Symbol("yu")
g2 = sympy.Symbol("g2")


def test_operators_expr():
    seesaw = v ** 2 / Λ
    assert seesaw == get_leading_mv(EFF_OPERATORS["1"])
    assert ye * loop * seesaw == get_leading_mv(EFF_OPERATORS["2"])
    assert yd * g2 * loop ** 2 * seesaw == get_leading_mv(EFF_OPERATORS["3a"])
    assert yd * loop * seesaw == get_leading_mv(EFF_OPERATORS["3b"])
    assert yu * loop * seesaw == get_leading_mv(EFF_OPERATORS["4a"])
    assert yu * g2 * loop ** 2 * seesaw == get_leading_mv(EFF_OPERATORS["4b"])
    assert yd * loop ** 2 * seesaw == get_leading_mv(EFF_OPERATORS["5a"])
    assert yd * loop * loopv2 * seesaw == get_leading_mv(EFF_OPERATORS["5b"])
    assert yd * g2 * loop ** 2 * loopv2 * seesaw == get_leading_mv(EFF_OPERATORS["5c"])
    assert yd * loop * loopv2 * seesaw == get_leading_mv(EFF_OPERATORS["5d"])
    assert yu * loop ** 2 * seesaw == get_leading_mv(EFF_OPERATORS["6a"])
    assert yu * loop * loopv2 * seesaw == get_leading_mv(EFF_OPERATORS["6b"])
    assert yu * loop * loopv2 * seesaw == get_leading_mv(EFF_OPERATORS["6c"])
    assert yu * g2 * loop ** 2 * loopv2 * seesaw == get_leading_mv(EFF_OPERATORS["6d"])
    assert ye * loop ** 2 * seesaw == get_leading_mv(EFF_OPERATORS["7"])
    assert ye * yu * yd * g2 * loop ** 2 * seesaw * loopv2 == get_leading_mv(
        EFF_OPERATORS["8"]
    )
    assert ye ** 2 * loop ** 2 * seesaw == get_leading_mv(EFF_OPERATORS["9"])
    assert ye * yd * loop ** 2 * seesaw == get_leading_mv(EFF_OPERATORS["10"])
    assert yd ** 2 * loop ** 3 * g2 * seesaw == get_leading_mv(EFF_OPERATORS["11a"])
    assert yd ** 2 * loop ** 2 * seesaw == get_leading_mv(EFF_OPERATORS["11b"])
    assert yu ** 2 * loop ** 2 * seesaw == get_leading_mv(EFF_OPERATORS["12a"])
    assert yu ** 2 * loop ** 3 * g2 * seesaw == get_leading_mv(EFF_OPERATORS["12b"])
    assert ye * yu * loop ** 2 * seesaw == get_leading_mv(EFF_OPERATORS["13"])
    assert yd * yu * loop ** 3 * g2 * seesaw == get_leading_mv(EFF_OPERATORS["14a"])
    assert yd * yu * loop ** 2 * seesaw == get_leading_mv(EFF_OPERATORS["14b"])
    assert yd * yu * loop ** 3 * g2 * seesaw == get_leading_mv(EFF_OPERATORS["15"])
    assert yd * yu * loop ** 3 * g2 * seesaw == get_leading_mv(EFF_OPERATORS["16"])
    assert yd * yu * loop ** 3 * g2 * seesaw == get_leading_mv(EFF_OPERATORS["17"])
    assert yd * yu * loop ** 3 * g2 * seesaw == get_leading_mv(EFF_OPERATORS["18"])
    assert yd ** 2 * yu * ye * g2 * loop ** 3 * loopv2 * seesaw == get_leading_mv(
        EFF_OPERATORS["19"]
    )
    assert yd * yu ** 2 * ye * g2 * loop ** 3 * loopv2 * seesaw == get_leading_mv(
        EFF_OPERATORS["20"]
    )
    yd ** 2 * ye ** 2 * yu ** 2 * g2 ** 2 * loop ** 4 * loopv2 ** 2 == get_leading_mv(
        EFF_OPERATORS["76"]
    )

    assert seesaw * loopv2 == get_leading_mv(EFF_OPERATORS["1p"])
    assert ye * loop * loopv2 * seesaw == get_leading_mv(EFF_OPERATORS["61a"])
    assert g2 * loop ** 2 * seesaw == get_leading_mv(DERIV_EFF_OPERATORS["D8i"])
    assert yd * yu * g2 * loop ** 2 * loopv2 * seesaw == get_leading_mv(
        DERIV_EFF_OPERATORS["D10b"]
    )


def test_operators_numerical():
    assert max(estimate_np_scale(EFF_OPERATORS["1"])) == 12
    assert max(estimate_np_scale(EFF_OPERATORS["2"])) == 8
    assert max(estimate_np_scale(EFF_OPERATORS["3a"])) == 4
    assert max(estimate_np_scale(EFF_OPERATORS["3b"])) == 8
    assert max(estimate_np_scale(EFF_OPERATORS["4a"])) == 10
    assert max(estimate_np_scale(EFF_OPERATORS["5a"])) == 6
    assert max(estimate_np_scale(EFF_OPERATORS["76"])) == -2
