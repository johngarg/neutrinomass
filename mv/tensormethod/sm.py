#!/usr/bin/env python3

from sympy import Rational
from core import Field

L = Field("L", dynkin="10001", charges={"y": Rational("-1/2")})
Q = Field("Q", dynkin="10101", charges={"y": Rational("1/6")})
H = Field("H", dynkin="00001", charges={"y": Rational("1/2")})
eb = Field("eb", dynkin="10000", charges={"y": Rational("1")})
ub = Field("ub", dynkin="10010", charges={"y": Rational("-2/3")})
db = Field("db", dynkin="10010", charges={"y": Rational("1/3")})
