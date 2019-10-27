#!/usr/bin/env python3

from sympy import Rational
from core import Field, FERMI, BOSE
import sympy.tensor.tensor as tensor

L = Field("L", dynkin="10001", charges={"y": Rational("-1/2")}, comm=FERMI)
Q = Field("Q", dynkin="10101", charges={"y": Rational("1/6")}, comm=FERMI)
H = Field("H", dynkin="00001", charges={"y": Rational("1/2")}, comm=BOSE)
eb = Field("eb", dynkin="10000", charges={"y": Rational("1")}, comm=FERMI)
ub = Field("ub", dynkin="10010", charges={"y": Rational("-2/3")}, comm=FERMI)
db = Field("db", dynkin="10010", charges={"y": Rational("1/3")}, comm=FERMI)

DL = Field("DL", dynkin="01001", charges={"y": Rational("-1/2")}, comm=FERMI)
DQ = Field("DQ", dynkin="01101", charges={"y": Rational("1/6")}, comm=FERMI)
DH = Field("DH", dynkin="11001", charges={"y": Rational("1/2")}, comm=BOSE)
Deb = Field("Deb", dynkin="01000", charges={"y": Rational("1")}, comm=FERMI)
Dub = Field("Dub", dynkin="01010", charges={"y": Rational("-2/3")}, comm=FERMI)
Ddb = Field("Ddb", dynkin="01010", charges={"y": Rational("1/3")}, comm=FERMI)
