#!/usr/bin/env python3

from sympy import Rational
import sympy.tensor.tensor as tensor

from neutrinomass.tensormethod.core import Field, FERMI, BOSE
from neutrinomass.tensormethod.lnv import BL_LIST

NF = 3

# fermions
L = Field("L", dynkin="10001", charges={"y": Rational("-1/2")}, comm=FERMI, nf=NF)
Q = Field("Q", dynkin="10101", charges={"y": Rational("1/6")}, comm=FERMI, nf=NF)
eb = Field("eb", dynkin="10000", charges={"y": Rational("1")}, comm=FERMI, nf=NF)
ub = Field("ub", dynkin="10010", charges={"y": Rational("-2/3")}, comm=FERMI, nf=NF)
db = Field("db", dynkin="10010", charges={"y": Rational("1/3")}, comm=FERMI, nf=NF)

# bosons
H = Field("H", dynkin="00001", charges={"y": Rational("1/2")}, comm=BOSE)
G = Field("G", dynkin="20110", comm=BOSE)
Gb = Field("Gb", dynkin="02110", comm=BOSE)
W = Field("W", dynkin="20002", comm=BOSE)
Wb = Field("Wb", dynkin="02002", comm=BOSE)
B = Field("B", dynkin="20000", comm=BOSE)
Bb = Field("Bb", dynkin="02000", comm=BOSE)

# set latex forms of base
L.latex = "L"
Q.latex = "Q"
H.latex = "H"
G.latex = "G"
W.latex = "W"
B.latex = "B"
eb.latex = r"\bar{e}"
ub.latex = r"\bar{u}"
db.latex = r"\bar{d}"
Gb.latex = r"\bar{G}"
Wb.latex = r"\bar{W}"
Bb.latex = r"\bar{B}"

# Add baryon number
L.charges["3b"] = 0
Q.charges["3b"] = 1
H.charges["3b"] = 0
eb.charges["3b"] = 0
ub.charges["3b"] = -1
db.charges["3b"] = -1

LNV_OPERATORS = {tuple([eval(field) for field in k]): v for k, v in BL_LIST.items()}
