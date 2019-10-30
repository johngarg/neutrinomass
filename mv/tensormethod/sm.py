#!/usr/bin/env python3

from sympy import Rational
from core import Field, FERMI, BOSE
from lnv import BL_LIST
import sympy.tensor.tensor as tensor

# fermions
L = Field("L", dynkin="10001", charges={"y": Rational("-1/2")}, comm=FERMI)
Q = Field("Q", dynkin="10101", charges={"y": Rational("1/6")}, comm=FERMI)
eb = Field("eb", dynkin="10000", charges={"y": Rational("1")}, comm=FERMI)
ub = Field("ub", dynkin="10010", charges={"y": Rational("-2/3")}, comm=FERMI)
db = Field("db", dynkin="10010", charges={"y": Rational("1/3")}, comm=FERMI)

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


def D(field, dynkin):
    """A derivative.

    Returns a new field with additional dotted and undotted indices.

    Example:
       >>> D(L, "01")
       DL(01001)(-1/2)

       >>> D(L, "21")
       DL(21001)(-1/2)

    """
    # derivative can only change one dotted and one undotted index
    assert abs(int(dynkin[0]) - field.dynkin_ints[0]) == 1
    assert abs(int(dynkin[1]) - field.dynkin_ints[1]) == 1

    # other info to construct field instance
    symbol = "D" + field.label
    new_field_dynkin = dynkin + field.dynkin[2:]
    rest = {"charges": field.charges, "comm": field.comm}

    new_field = Field(symbol, dynkin=new_field_dynkin, **rest)
    return new_field


LNV_OPERATORS = {tuple([eval(field) for field in k]): v for k, v in BL_LIST.items()}
