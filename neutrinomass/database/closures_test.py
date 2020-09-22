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

SEESAW = v ** 2 / Λ

dGJ_RESULTS = {
    "1": SEESAW,
    "2": SEESAW * ye * loop,
    "3a": SEESAW * yd * loop ** 2 * g2,
    "3b": SEESAW * yd * loop,
    "4a": SEESAW * yu * loop,
    "4b": SEESAW * yu * loop ** 2 * g2,
    "5a": SEESAW * yd * loop ** 2,
    "6a": SEESAW * yu * loop ** 2,
    "7": SEESAW * ye * loop ** 2 * g2 * loopv2,
    "8": SEESAW * ye * yd * yu * loop ** 2,
    "9": SEESAW * ye ** 2 * loop ** 2,
    "10": SEESAW * ye * yd * loop ** 2,
    "11a": SEESAW * yd ** 2 * g2 * loop ** 3,
    "11b": SEESAW * yd ** 2 * loop ** 2,
    "12b": SEESAW * yu ** 2 * loop ** 2,
    "12b": SEESAW * yu ** 2 * g2 * loop ** 3,
    "13": SEESAW * ye * yu * loop ** 2,
    "14a": SEESAW * yd * yu * g2 * loop ** 3,
    "14b": SEESAW * yd * yu * loop ** 2,
    "15": SEESAW * yd * yu * g2 * loop ** 3,
    "16": SEESAW * yd * yu * g2 ** 2 * loop ** 4,
    "17": SEESAW * yd * yu * g2 ** 2 * loop ** 4,
    "18": SEESAW * yd * yu * g2 ** 2 * loop ** 4,
    "19": SEESAW * yd ** 2 * yu * ye * loop ** 3,
    "20": SEESAW * yd * yu ** 2 * ye * loop ** 3,
    "21a": SEESAW * ye * yu * loop ** 2 * loopv2,
    "21b": SEESAW * ye * yu * loop ** 2 * loopv2,
    "22a": SEESAW * g2 * loop ** 3,
    "23a": SEESAW * loopv2 * ye * yd * loop ** 2,
    "24a": SEESAW * yd ** 2 * loop ** 3,
    "24b": SEESAW * yd ** 2 * loop ** 3,
    "25a": SEESAW * yd * yu * loop ** 2 * loopv2,
    "26a": SEESAW * yd * ye * loop ** 3,
    "26b": SEESAW * yd * ye * loop ** 2 * loopv2,
    "27a": SEESAW * loop ** 3 * g2,
    "27b": SEESAW * loop ** 3 * g2,
    "28a": SEESAW * yd * yu * loop ** 3,
    "28b": SEESAW * yd * yu * loop ** 3,
    "28c": SEESAW * yd * yu * loop ** 3,
    "29a": SEESAW * yu ** 2 * loop ** 2 * loopv2,
    "29b": SEESAW * g2 * loop ** 3,
    "30a": SEESAW * ye * yu * loop ** 3,
    "30b": SEESAW * ye * yu * loop ** 2 * loopv2,
    "31a": SEESAW * yd * yu * loop ** 2 * loopv2,
    "31b": SEESAW * yd * yu * loop ** 2 * loopv2,
    "37": SEESAW * ye ** 2 * loop ** 5 * yd ** 2 * g2,
    "43a": SEESAW * g2 * loop ** 4 * yu * yd,
    "43b": SEESAW * g2 * loop ** 4 * yu * yd,
    "43c": SEESAW * g2 * loop ** 4 * yu * yd,
    "49": g2 * loop ** 3 * SEESAW,
    "50a": yd * yu * g2 * loop ** 3 * SEESAW,
    "52a": yd * yu * g2 * loop ** 4 * SEESAW,
    "53": yd ** 2 * yu ** 2 * g2 * loop ** 5 * SEESAW,
    "57": ye * yu * g2 * loop ** 4 * SEESAW,
    "59a": ye * yu * yd ** 2 * loop ** 4 * SEESAW,
    "60a": ye * yu ** 2 * yd * loop ** 4 * SEESAW,
    "65a": yu * yd * loop ** 4 * g2 * SEESAW,
    "75": ye * yd * yu ** 2 * loop ** 3 * loopv2 * SEESAW,
}


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
    assert dGJ_RESULTS["21a"] == get_leading_mv(EFF_OPERATORS["21a"])
    assert dGJ_RESULTS["21b"] == get_leading_mv(EFF_OPERATORS["21b"])
    assert loop ** 4 * yd * ye ** 2 * yu * seesaw == get_leading_mv(EFF_OPERATORS["37"])
    assert g2 * loop ** 3 * loopv2 * yd * yu * seesaw == get_leading_mv(
        EFF_OPERATORS["43a"]
    )
    assert g2 * loop ** 3 * loopv2 * yd * yu * seesaw == get_leading_mv(
        EFF_OPERATORS["50a"]
    )
    assert g2 * loop ** 3 * loopv2 * yd * yu * seesaw == get_leading_mv(
        EFF_OPERATORS["52a"]
    )
    assert (
        g2 ** 2 * loop ** 4 * loopv2 ** 2 * yd ** 2 * yu ** 2 * seesaw
        == get_leading_mv(EFF_OPERATORS["53"])
    )
    assert (
        seesaw * yd ** 2 * ye ** 2 * yu ** 2 * g2 ** 2 * loop ** 4 * loopv2 ** 2
        == get_leading_mv(EFF_OPERATORS["76"])
    )
    assert loop ** 3 * ye * yu * seesaw == get_leading_mv(EFF_OPERATORS["57"])
    assert g2 * loop ** 3 * loopv2 ** 2 * ye * yu * yd ** 2 * seesaw == get_leading_mv(
        EFF_OPERATORS["59a"]
    )
    assert g2 * loop ** 4 * loopv2 * ye * yu ** 2 * yd * seesaw == get_leading_mv(
        EFF_OPERATORS["60a"]
    )
    assert loop ** 3 * yd * ye * seesaw == get_leading_mv(EFF_OPERATORS["75"])

    assert seesaw * loopv2 == get_leading_mv(EFF_OPERATORS["1p"])
    assert ye * loop * loopv2 * seesaw == get_leading_mv(EFF_OPERATORS["61a"])
    assert g2 * loop ** 2 * seesaw == get_leading_mv(DERIV_EFF_OPERATORS["D8i"])
    assert yd * yu * g2 * loop ** 2 * loopv2 * seesaw == get_leading_mv(
        DERIV_EFF_OPERATORS["D10b"]
    )


def test_operators_numerical():
    from math import log10

    proc = lambda n: round(log10(max(n)))

    assert proc(estimate_np_scale(EFF_OPERATORS["1"])) == 12
    assert proc(estimate_np_scale(EFF_OPERATORS["2"])) == 8
    assert proc(estimate_np_scale(EFF_OPERATORS["3a"])) == 5
    assert proc(estimate_np_scale(EFF_OPERATORS["3b"])) == 8
    assert proc(estimate_np_scale(EFF_OPERATORS["4a"])) == 10
    assert proc(estimate_np_scale(EFF_OPERATORS["5a"])) == 6
    assert proc(estimate_np_scale(EFF_OPERATORS["76"])) == -2


def test_dgj():
    checked = {
        # vanishing loop
        "7",
        "16",
        "17",
        "18",
        "22a",
        "27a",
        "27b",
        "29a",
        "29b",
        "49",
        "50a",
        "52a",
        "57",
        "75",
        # vanishing Higgs combo
        "8",
        "19",
        "20",
        "28a",
        "28b",
        "28c",
        "43a",
        "43b",
        "43c",
        "53",
        "59a",
        "60a",
        "65a",
        # unclear
        "37",
    }
    for k, v in dGJ_RESULTS.items():
        if k in checked:
            continue

        expr = get_leading_mv(EFF_OPERATORS[k])
        assert v == expr
        # if v != expr:
        #     print(f"{k}: {expr}")
