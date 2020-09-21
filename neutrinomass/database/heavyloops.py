#!/usr/bin/env python3

"""Tables of the fields in arXiv:1204.5862 `Systematic study of the $d=5$
Weinberg operator at one-loop order' for filtering

"""

from sympy import Rational, Symbol
from functools import reduce
import math

# su2
tbl_2_su2_reps_func = lambda m: [
    [m, m + 1, m, m + 1],
    [m + 1, m, m + 1, m],
    [m, m + 1, m + 2, m + 1],
    [m + 1, m, m + 1, m + 2],
    [m + 1, m + 2, m + 1, m],
    [m + 2, m + 1, m, m + 1],
    [m + 1, m + 2, m + 1, m + 2],
    [m + 2, m + 1, m + 2, m + 1],
]

tbl_3_su2_reps_func = lambda m: [
    [m, m + 2, m + 1],
    [m + 1, m + 1, m],
    [m + 1, m + 1, m + 2],
    [m + 2, m, m + 1],
    [m + 2, m + 2, m + 1],
]

tbl_2_su2_reps = (
    tbl_2_su2_reps_func(1)
    + tbl_2_su2_reps_func(2)
    + tbl_2_su2_reps_func(3)
    + tbl_2_su2_reps_func(4)
    + tbl_2_su2_reps_func(5)
)
tbl_3_su2_reps = (
    tbl_3_su2_reps_func(1)
    + tbl_3_su2_reps_func(2)
    + tbl_3_su2_reps_func(3)
    + tbl_3_su2_reps_func(4)
    + tbl_3_su2_reps_func(5)
)


def tbl_2_1(i: list, a: Rational):
    return [
        (
            f"S,00,{i[0]-1},{a}",
            f"S,00,{i[1]-1},{a-Rational(1,2)}",
            f"F,00,{i[2]-1},{a}",
            f"S,00,{i[3]-1},{Rational(1,2) + a}",
        ),
        (
            f"S,10,{i[0]-1},{a}",
            f"S,01,{i[1]-1},{a-Rational(1,2)}",
            f"F,10,{i[2]-1},{a}",
            f"S,01,{i[3]-1},{Rational(1,2) + a}",
        ),
        (
            f"S,01,{i[0]-1},{a}",
            f"S,10,{i[1]-1},{a-Rational(1,2)}",
            f"F,01,{i[2]-1},{a}",
            f"S,10,{i[3]-1},{Rational(1,2) + a}",
        ),
        (
            f"S,02,{i[0]-1},{a}",
            f"S,20,{i[1]-1},{a-Rational(1,2)}",
            f"F,02,{i[2]-1},{a}",
            f"S,20,{i[3]-1},{Rational(1,2) + a}",
        ),
        (
            f"S,20,{i[0]-1},{a}",
            f"S,02,{i[1]-1},{a-Rational(1,2)}",
            f"F,20,{i[2]-1},{a}",
            f"S,02,{i[3]-1},{Rational(1,2) + a}",
        ),
        (
            f"S,11,{i[0]-1},{a}",
            f"S,11,{i[1]-1},{a-Rational(1,2)}",
            f"F,11,{i[2]-1},{a}",
            f"S,11,{i[3]-1},{Rational(1,2) + a}",
        ),
    ]


def tbl_2_2(i: list, a: Rational):
    return [
        (
            f"F,00,{i[0]-1},{a}",
            f"S,00,{i[1]-1},{a+Rational(1,2)}",
            f"S,00,{i[2]-1},{a}",
            f"F,00,{i[3]-1},{Rational(1,2) + a}",
        ),
        (
            f"F,10,{i[0]-1},{a}",
            f"S,01,{i[1]-1},{a+Rational(1,2)}",
            f"S,10,{i[2]-1},{a}",
            f"F,01,{i[3]-1},{Rational(1,2) + a}",
        ),
        (
            f"F,01,{i[0]-1},{a}",
            f"S,10,{i[1]-1},{a+Rational(1,2)}",
            f"S,01,{i[2]-1},{a}",
            f"F,10,{i[3]-1},{Rational(1,2) + a}",
        ),
        (
            f"F,02,{i[0]-1},{a}",
            f"S,20,{i[1]-1},{a+Rational(1,2)}",
            f"S,02,{i[2]-1},{a}",
            f"F,20,{i[3]-1},{Rational(1,2) + a}",
        ),
        (
            f"F,20,{i[0]-1},{a}",
            f"S,02,{i[1]-1},{a+Rational(1,2)}",
            f"S,20,{i[2]-1},{a}",
            f"F,02,{i[3]-1},{Rational(1,2) + a}",
        ),
        (
            f"F,11,{i[0]-1},{a}",
            f"S,11,{i[1]-1},{a+Rational(1,2)}",
            f"S,11,{i[2]-1},{a}",
            f"F,11,{i[3]-1},{Rational(1,2) + a}",
        ),
    ]


def tbl_2_3(i: list, a: Rational):
    return [
        (
            f"F,00,{i[0]-1},{a}",
            f"F,00,{i[1]-1},{a+Rational(1,2)}",
            f"S,00,{i[2]-1},{a}",
            f"F,00,{i[3]-1},{-Rational(1,2) + a}",
        ),
        (
            f"F,10,{i[0]-1},{a}",
            f"F,01,{i[1]-1},{a+Rational(1,2)}",
            f"S,10,{i[2]-1},{a}",
            f"F,01,{i[3]-1},{-Rational(1,2) + a}",
        ),
        (
            f"F,01,{i[0]-1},{a}",
            f"F,10,{i[1]-1},{a+Rational(1,2)}",
            f"S,01,{i[2]-1},{a}",
            f"F,10,{i[3]-1},{-Rational(1,2) + a}",
        ),
        (
            f"F,02,{i[0]-1},{a}",
            f"F,20,{i[1]-1},{a+Rational(1,2)}",
            f"S,02,{i[2]-1},{a}",
            f"F,20,{i[3]-1},{-Rational(1,2) + a}",
        ),
        (
            f"F,20,{i[0]-1},{a}",
            f"F,02,{i[1]-1},{a+Rational(1,2)}",
            f"S,20,{i[2]-1},{a}",
            f"F,02,{i[3]-1},{-Rational(1,2) + a}",
        ),
        (
            f"F,11,{i[0]-1},{a}",
            f"F,11,{i[1]-1},{a+Rational(1,2)}",
            f"S,11,{i[2]-1},{a}",
            f"F,11,{i[3]-1},{-Rational(1,2) + a}",
        ),
    ]


def tbl_3_1(i: list, a: Rational):
    return [
        (
            f"S,00,{i[0]-1},{a}",
            f"S,00,{i[1]-1},{1+a}",
            f"F,00,{i[2]-1},{Rational(1,2) + a}",
        ),
        (
            f"S,10,{i[0]-1},{a}",
            f"S,10,{i[1]-1},{1+a}",
            f"F,01,{i[2]-1},{Rational(1,2) + a}",
        ),
        (
            f"S,01,{i[0]-1},{a}",
            f"S,01,{i[1]-1},{1+a}",
            f"F,10,{i[2]-1},{Rational(1,2) + a}",
        ),
        (
            f"S,02,{i[0]-1},{a}",
            f"S,02,{i[1]-1},{1+a}",
            f"F,20,{i[2]-1},{Rational(1,2) + a}",
        ),
        (
            f"S,20,{i[0]-1},{a}",
            f"S,20,{i[1]-1},{1+a}",
            f"F,02,{i[2]-1},{Rational(1,2) + a}",
        ),
        (
            f"S,11,{i[0]-1},{a}",
            f"S,11,{i[1]-1},{1+a}",
            f"F,11,{i[2]-1},{Rational(1,2) + a}",
        ),
    ]


def generate_models():
    # Maximum hypercharge in models is 3
    counter = -3 - Rational("1/6")
    y_range = []
    for i in range(37):
        counter += Rational("1/6")
        y_range.append(counter)

    models = []
    for y in y_range:
        models += reduce(
            lambda x, y: x + y, [tbl_2_1(struct, y) for struct in tbl_2_su2_reps]
        )
        models += reduce(
            lambda x, y: x + y, [tbl_2_2(struct, y) for struct in tbl_2_su2_reps]
        )
        models += reduce(
            lambda x, y: x + y, [tbl_2_3(struct, y) for struct in tbl_2_su2_reps]
        )
        models += reduce(
            lambda x, y: x + y, [tbl_3_1(struct, y) for struct in tbl_3_su2_reps]
        )

    return models
