#!/usr/bin/env python3

"""Results from the Hilbert Series for the SMEFT up to dimension 7 and for the
ΔL = 2 SMEFT up to dimension 11. Derivatives kept up until dimension 9."""

from sympy import Function
from sympy import symbols
from sympy.abc import X

# Nf = number of fermion generations
# D = symbolic representation for the derivative
Nf, D = symbols("Nf D")

# Treat fields as sympy functions to easily deal with derivatives
H, Hd, L, Ld, Q, Qd = symbols("H Hd L Ld Q Qd", cls=Function)
eb, ebd, ub, ubd, db, dbd = symbols("eb ebd ub ubd db dbd", cls=Function)
G, Gb, W, Wb, B, Bb = symbols("G Gb W Wb B Bb", cls=Function)

# Dimension 6 HS with field strengths and IBP + EOM redundancies removed
# H6 = (
#     (Nf ** 2 * db(X) ** 2 * dbd(X) ** 2) / 2
#     + (Nf ** 4 * db(X) ** 2 * dbd(X) ** 2) / 2
#     + Nf ** 4 * db(X) * dbd(X) * eb(X) * ebd(X)
#     + -(Nf ** 2 * eb(X) ** 2 * ebd(X) ** 2) / 4
#     + (Nf ** 3 * eb(X) ** 2 * ebd(X) ** 2) / 2
#     + (Nf ** 4 * eb(X) ** 2 * ebd(X) ** 2) / 4
#     + -G(X) ** 3
#     + Gb(X) ** 3
#     + B(X) ** 2 * H(X) * Hd(X)
#     + Bb(X) ** 2 * H(X) * Hd(X)
#     + D * Nf ** 2 * db(X) * dbd(X) * H(X) * Hd(X)
#     + -D * Nf ** 2 * eb(X) * ebd(X) * H(X) * Hd(X)
#     + G(X) ** 2 * H(X) * Hd(X)
#     + Gb(X) ** 2 * H(X) * Hd(X)
#     + 2 * D ** 2 * H(X) ** 2 * Hd(X) ** 2
#     + -H(X) ** 3 * Hd(X) ** 3
#     + Nf ** 2 * B(X) * eb(X) * Hd(X) * L(X)
#     + Nf ** 2 * eb(X) * H(X) * Hd(X) ** 2 * L(X)
#     + -Nf ** 2 * Bb(X) * ebd(X) * H(X) * Ld(X)
#     + Nf ** 2 * ebd(X) * H(X) ** 2 * Hd(X) * Ld(X)
#     + Nf ** 4 * db(X) * dbd(X) * L(X) * Ld(X)
#     + -Nf ** 4 * eb(X) * ebd(X) * L(X) * Ld(X)
#     + 2 * D * Nf ** 2 * H(X) * Hd(X) * L(X) * Ld(X)
#     + (Nf ** 2 * L(X) ** 2 * Ld(X) ** 2) / 2
#     + -(Nf ** 4 * L(X) ** 2 * Ld(X) ** 2) / 2
#     + Nf ** 2 * B(X) * db(X) * Hd(X) * Q(X)
#     + Nf ** 2 * db(X) * G(X) * Hd(X) * Q(X)
#     + -Nf ** 2 * db(X) * H(X) * Hd(X) ** 2 * Q(X)
#     + Nf ** 4 * db(X) * ebd(X) * Ld(X) * Q(X)
#     + (Nf ** 2 * L(X) * Q(X) ** 3) / 3
#     + -(2 * Nf ** 4 * L(X) * Q(X) ** 3) / 3
#     + Nf ** 2 * Bb(X) * dbd(X) * H(X) * Qd(X)
#     + Nf ** 2 * dbd(X) * Gb(X) * H(X) * Qd(X)
#     + -Nf ** 2 * dbd(X) * H(X) ** 2 * Hd(X) * Qd(X)
#     + Nf ** 4 * dbd(X) * eb(X) * L(X) * Qd(X)
#     + 2 * Nf ** 4 * db(X) * dbd(X) * Q(X) * Qd(X)
#     + -Nf ** 4 * eb(X) * ebd(X) * Q(X) * Qd(X)
#     + 2 * D * Nf ** 2 * H(X) * Hd(X) * Q(X) * Qd(X)
#     + 2 * Nf ** 4 * L(X) * Ld(X) * Q(X) * Qd(X)
#     + -Nf ** 2 * Q(X) ** 2 * Qd(X) ** 2
#     + Nf ** 4 * Q(X) ** 2 * Qd(X) ** 2
#     + (Nf ** 2 * Ld(X) * Qd(X) ** 3) / 3
#     + -(2 * Nf ** 4 * Ld(X) * Qd(X) ** 3) / 3
#     + D * Nf ** 2 * dbd(X) * H(X) ** 2 * ub(X)
#     + Nf ** 2 * B(X) * H(X) * Q(X) * ub(X)
#     + -Nf ** 2 * G(X) * H(X) * Q(X) * ub(X)
#     + Nf ** 2 * H(X) ** 2 * Hd(X) * Q(X) * ub(X)
#     + 2 * Nf ** 4 * eb(X) * L(X) * Q(X) * ub(X)
#     + -2 * Nf ** 4 * db(X) * Q(X) ** 2 * ub(X)
#     + Nf ** 4 * db(X) * Ld(X) * Qd(X) * ub(X)
#     + (Nf ** 3 * eb(X) * Qd(X) ** 2 * ub(X)) / 2
#     + -(Nf ** 4 * eb(X) * Qd(X) ** 2 * ub(X)) / 2
#     + Nf ** 4 * db(X) * eb(X) * ub(X) ** 2
#     + D * Nf ** 2 * db(X) * Hd(X) ** 2 * ubd(X)
#     + -Nf ** 4 * dbd(X) * L(X) * Q(X) * ubd(X)
#     + (Nf ** 3 * ebd(X) * Q(X) ** 2 * ubd(X)) / 2
#     + (Nf ** 4 * ebd(X) * Q(X) ** 2 * ubd(X)) / 2
#     + -Nf ** 2 * Bb(X) * Hd(X) * Qd(X) * ubd(X)
#     + Nf ** 2 * Gb(X) * Hd(X) * Qd(X) * ubd(X)
#     + Nf ** 2 * H(X) * Hd(X) ** 2 * Qd(X) * ubd(X)
#     + -2 * Nf ** 4 * ebd(X) * Ld(X) * Qd(X) * ubd(X)
#     + 2 * Nf ** 4 * dbd(X) * Qd(X) ** 2 * ubd(X)
#     + 2 * Nf ** 4 * db(X) * dbd(X) * ub(X) * ubd(X)
#     + -Nf ** 4 * eb(X) * ebd(X) * ub(X) * ubd(X)
#     + D * Nf ** 2 * H(X) * Hd(X) * ub(X) * ubd(X)
#     + Nf ** 4 * L(X) * Ld(X) * ub(X) * ubd(X)
#     + -2 * Nf ** 4 * Q(X) * Qd(X) * ub(X) * ubd(X)
#     + Nf ** 4 * dbd(X) * ebd(X) * ubd(X) ** 2
#     + (Nf ** 2 * ub(X) ** 2 * ubd(X) ** 2) / 2
#     + -(Nf ** 4 * ub(X) ** 2 * ubd(X) ** 2) / 2
#     + B(X) * H(X) * Hd(X) * W(X)
#     + Nf ** 2 * eb(X) * Hd(X) * L(X) * W(X)
#     + -Nf ** 2 * db(X) * Hd(X) * Q(X) * W(X)
#     + Nf ** 2 * H(X) * Q(X) * ub(X) * W(X)
#     + H(X) * Hd(X) * W(X) ** 2
#     + W(X) ** 3
#     + -Bb(X) * H(X) * Hd(X) * Wb(X)
#     + Nf ** 2 * ebd(X) * H(X) * Ld(X) * Wb(X)
#     + Nf ** 2 * dbd(X) * H(X) * Qd(X) * Wb(X)
#     + -Nf ** 2 * Hd(X) * Qd(X) * ubd(X) * Wb(X)
#     + H(X) * Hd(X) * Wb(X) ** 2
#     + Wb(X) ** 3
# )

H6 = (
    db(X) ** 2 * dbd(X) ** 2,
    db(X) * dbd(X) * eb(X) * ebd(X),
    eb(X) ** 2 * ebd(X) ** 2,
    H(X) ** 3 * Hd(X) ** 3,
    -eb(X) * H(X) * Hd(X) ** 2 * L(X),
    ebd(X) * H(X) ** 2 * Hd(X) * Ld(X),
    db(X) * dbd(X) * L(X) * Ld(X),
    eb(X) * ebd(X) * L(X) * Ld(X),
    -L(X) ** 2 * Ld(X) ** 2,
    db(X) * H(X) * Hd(X) ** 2 * Q(X),
    db(X) * ebd(X) * Ld(X) * Q(X),
    L(X) * Q(X) ** 3,
    dbd(X) * H(X) ** 2 * Hd(X) * Qd(X),
    -dbd(X) * eb(X) * L(X) * Qd(X),
    2 * db(X) * dbd(X) * Q(X) * Qd(X),
    eb(X) * ebd(X) * Q(X) * Qd(X),
    2 * L(X) * Ld(X) * Q(X) * Qd(X),
    -2 * Q(X) ** 2 * Qd(X) ** 2,
    Ld(X) * Qd(X) ** 3,
    3 * D * db(X) * dbd(X) * H(X) * Hd(X),
    3 * D * eb(X) * ebd(X) * H(X) * Hd(X),
    -6 * D * H(X) * Hd(X) * L(X) * Ld(X),
    6 * D * H(X) * Hd(X) * Q(X) * Qd(X),
    4 * D ** 2 * H(X) ** 2 * Hd(X) ** 2,
    4 * D ** 2 * eb(X) * Hd(X) * L(X),
    -4 * D ** 2 * ebd(X) * H(X) * Ld(X),
    4 * D ** 2 * db(X) * Hd(X) * Q(X),
    4 * D ** 2 * dbd(X) * H(X) * Qd(X),
    D ** 3 * db(X) * dbd(X),
    -D ** 3 * eb(X) * ebd(X),
    D ** 3 * L(X) * Ld(X),
    D ** 3 * Q(X) * Qd(X),
    D ** 4 * H(X) * Hd(X),
    H(X) ** 2 * Hd(X) * Q(X) * ub(X),
    -2 * eb(X) * L(X) * Q(X) * ub(X),
    2 * db(X) * Q(X) ** 2 * ub(X),
    db(X) * Ld(X) * Qd(X) * ub(X),
    eb(X) * Qd(X) ** 2 * ub(X),
    -D * dbd(X) * H(X) ** 2 * ub(X),
    4 * D ** 2 * H(X) * Q(X) * ub(X),
    db(X) * eb(X) * ub(X) ** 2,
    dbd(X) * L(X) * Q(X) * ubd(X),
    -ebd(X) * Q(X) ** 2 * ubd(X),
    H(X) * Hd(X) ** 2 * Qd(X) * ubd(X),
    2 * ebd(X) * Ld(X) * Qd(X) * ubd(X),
    2 * dbd(X) * Qd(X) ** 2 * ubd(X),
    -D * db(X) * Hd(X) ** 2 * ubd(X),
    4 * D ** 2 * Hd(X) * Qd(X) * ubd(X),
    2 * db(X) * dbd(X) * ub(X) * ubd(X),
    eb(X) * ebd(X) * ub(X) * ubd(X),
    -L(X) * Ld(X) * ub(X) * ubd(X),
    2 * Q(X) * Qd(X) * ub(X) * ubd(X),
    3 * D * H(X) * Hd(X) * ub(X) * ubd(X),
    D ** 3 * ub(X) * ubd(X),
    -dbd(X) * ebd(X) * ubd(X) ** 2,
    ub(X) ** 2 * ubd(X) ** 2,
)

H7_LNV = (
    Nf ** 2 * D * ebd(X) * H(X) ** 3 * L(X)
    + Nf * D ** 2 * H(X) ** 2 * L(X) ** 2
    + Nf ** 2 * D ** 2 * H(X) ** 2 * L(X) ** 2
    - (Nf * B(X) * H(X) ** 2 * L(X) ** 2) / 2
    + -(Nf ** 2 * B(X) * H(X) ** 2 * L(X) ** 2) / 2
    + (Nf * H(X) ** 3 * Hd(X) * L(X) ** 2) / 2
    + (Nf ** 2 * H(X) ** 3 * Hd(X) * L(X) ** 2) / 2
    + -(Nf ** 2 * eb(X) * H(X) * L(X) ** 3) / 3
    + (2 * Nf ** 4 * eb(X) * H(X) * L(X) ** 3) / 3
    + 2 * Nf ** 4 * db(X) * H(X) * L(X) ** 2 * Q(X)
    + -Nf ** 4 * db(X) * ebd(X) * H(X) * L(X) * ubd(X)
    + (Nf ** 3 * D * db(X) * L(X) ** 2 * ubd(X)) / 2
    + (Nf ** 4 * D * db(X) * L(X) ** 2 * ubd(X)) / 2
    + -Nf ** 4 * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + Nf ** 2 * H(X) ** 2 * L(X) ** 2 * W(X)
)

H9_LNV = (
    (Nf * D ** 2 * ebd(X) ** 2 * H(X) ** 4) / 2
    + (Nf ** 2 * D ** 2 * ebd(X) ** 2 * H(X) ** 4) / 2
    + 3 * Nf ** 2 * D ** 3 * ebd(X) * H(X) ** 3 * L(X)
    + -Nf ** 2 * D * B(X) * ebd(X) * H(X) ** 3 * L(X)
    + Nf ** 2 * D * Bb(X) * ebd(X) * H(X) ** 3 * L(X)
    + Nf ** 2 * D * ebd(X) * H(X) ** 4 * Hd(X) * L(X)
    + -(3 * Nf * D ** 4 * H(X) ** 2 * L(X) ** 2) / 2
    + (3 * Nf ** 2 * D ** 4 * H(X) ** 2 * L(X) ** 2) / 2
    - (Nf * D ** 2 * B(X) * H(X) ** 2 * L(X) ** 2) / 2
    + -(7 * Nf ** 2 * D ** 2 * B(X) * H(X) ** 2 * L(X) ** 2) / 2
    + (Nf * B(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    + (Nf ** 2 * B(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    - -Nf * D ** 2 * Bb(X) * H(X) ** 2 * L(X) ** 2
    + 2 * Nf ** 2 * D ** 2 * Bb(X) * H(X) ** 2 * L(X) ** 2
    + (Nf * Bb(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    + -(Nf ** 2 * Bb(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    - (Nf ** 3 * D * db(X) * dbd(X) * H(X) ** 2 * L(X) ** 2) / 2
    + -(5 * Nf ** 4 * D * db(X) * dbd(X) * H(X) ** 2 * L(X) ** 2) / 2
    - (Nf ** 3 * D * eb(X) * ebd(X) * H(X) ** 2 * L(X) ** 2) / 2
    + -(5 * Nf ** 4 * D * eb(X) * ebd(X) * H(X) ** 2 * L(X) ** 2) / 2
    + (Nf * G(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    + (Nf ** 2 * G(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    + -(Nf * Gb(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    + (Nf ** 2 * Gb(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    + (Nf * D ** 2 * H(X) ** 3 * Hd(X) * L(X) ** 2) / 2
    + -(11 * Nf ** 2 * D ** 2 * H(X) ** 3 * Hd(X) * L(X) ** 2) / 2
    - (Nf * B(X) * H(X) ** 3 * Hd(X) * L(X) ** 2) / 2
    + -(Nf ** 2 * B(X) * H(X) ** 3 * Hd(X) * L(X) ** 2) / 2
    + (Nf * H(X) ** 4 * Hd(X) ** 2 * L(X) ** 2) / 2
    + (Nf ** 2 * H(X) ** 4 * Hd(X) ** 2 * L(X) ** 2) / 2
    - -(Nf ** 2 * D ** 2 * eb(X) * H(X) * L(X) ** 3) / 3
    + (10 * Nf ** 4 * D ** 2 * eb(X) * H(X) * L(X) ** 3) / 3
    + Nf ** 4 * B(X) * eb(X) * H(X) * L(X) ** 3
    + -Nf ** 4 * eb(X) * H(X) ** 2 * Hd(X) * L(X) ** 3
    + (13 * Nf ** 3 * eb(X) ** 2 * L(X) ** 4) / 24
    + (7 * Nf ** 4 * eb(X) ** 2 * L(X) ** 4) / 24
    - -(Nf ** 5 * eb(X) ** 2 * L(X) ** 4) / 24
    + (5 * Nf ** 6 * eb(X) ** 2 * L(X) ** 4) / 24
    + (Nf ** 3 * ebd(X) * H(X) ** 3 * L(X) ** 2 * Ld(X)) / 2
    + -(Nf ** 4 * ebd(X) * H(X) ** 3 * L(X) ** 2 * Ld(X)) / 2
    + (Nf ** 2 * D * H(X) ** 2 * L(X) ** 3 * Ld(X)) / 3
    - (Nf ** 3 * D * H(X) ** 2 * L(X) ** 3 * Ld(X)) / 2
    + -(13 * Nf ** 4 * D * H(X) ** 2 * L(X) ** 3 * Ld(X)) / 6
    + 5 * Nf ** 4 * D * db(X) * ebd(X) * H(X) ** 2 * L(X) * Q(X)
    + -10 * Nf ** 4 * D ** 2 * db(X) * H(X) * L(X) ** 2 * Q(X)
    + 3 * Nf ** 4 * B(X) * db(X) * H(X) * L(X) ** 2 * Q(X)
    + 3 * Nf ** 4 * db(X) * G(X) * H(X) * L(X) ** 2 * Q(X)
    + -3 * Nf ** 4 * db(X) * H(X) ** 2 * Hd(X) * L(X) ** 2 * Q(X)
    + (Nf ** 4 * db(X) * eb(X) * L(X) ** 3 * Q(X)) / 3
    + -(5 * Nf ** 6 * db(X) * eb(X) * L(X) ** 3 * Q(X)) / 3
    + (3 * Nf ** 3 * db(X) ** 2 * L(X) ** 2 * Q(X) ** 2) / 2
    + -(5 * Nf ** 6 * db(X) ** 2 * L(X) ** 2 * Q(X) ** 2) / 2
    + (Nf ** 3 * dbd(X) * H(X) ** 3 * L(X) ** 2 * Qd(X)) / 2
    + -(Nf ** 4 * dbd(X) * H(X) ** 3 * L(X) ** 2 * Qd(X)) / 2
    + Nf ** 4 * ebd(X) * H(X) ** 3 * L(X) * Q(X) * Qd(X)
    - -(Nf ** 3 * D * H(X) ** 2 * L(X) ** 2 * Q(X) * Qd(X)) / 2
    + (13 * Nf ** 4 * D * H(X) ** 2 * L(X) ** 2 * Q(X) * Qd(X)) / 2
    + -Nf ** 4 * H(X) ** 3 * L(X) ** 2 * Q(X) * ub(X)
    + Nf ** 4 * D * db(X) * ebd(X) ** 2 * H(X) ** 2 * ubd(X)
    + -7 * Nf ** 4 * D ** 2 * db(X) * ebd(X) * H(X) * L(X) * ubd(X)
    + Nf ** 4 * B(X) * db(X) * ebd(X) * H(X) * L(X) * ubd(X)
    + -Nf ** 4 * Bb(X) * db(X) * ebd(X) * H(X) * L(X) * ubd(X)
    + Nf ** 4 * db(X) * ebd(X) * G(X) * H(X) * L(X) * ubd(X)
    + -Nf ** 4 * db(X) * ebd(X) * Gb(X) * H(X) * L(X) * ubd(X)
    + Nf ** 4 * db(X) * ebd(X) * H(X) ** 2 * Hd(X) * L(X) * ubd(X)
    + -Nf ** 4 * D ** 3 * db(X) * L(X) ** 2 * ubd(X)
    + 2 * Nf ** 4 * D * B(X) * db(X) * L(X) ** 2 * ubd(X)
    + (Nf ** 3 * D * Bb(X) * db(X) * L(X) ** 2 * ubd(X)) / 2
    + -(3 * Nf ** 4 * D * Bb(X) * db(X) * L(X) ** 2 * ubd(X)) / 2
    + Nf ** 6 * db(X) ** 2 * dbd(X) * L(X) ** 2 * ubd(X)
    + -Nf ** 6 * db(X) * eb(X) * ebd(X) * L(X) ** 2 * ubd(X)
    + 2 * Nf ** 4 * D * db(X) * G(X) * L(X) ** 2 * ubd(X)
    + -(Nf ** 3 * D * db(X) * Gb(X) * L(X) ** 2 * ubd(X)) / 2
    + (3 * Nf ** 4 * D * db(X) * Gb(X) * L(X) ** 2 * ubd(X)) / 2
    + -5 * Nf ** 4 * D * db(X) * H(X) * Hd(X) * L(X) ** 2 * ubd(X)
    + (Nf ** 4 * db(X) * L(X) ** 3 * Ld(X) * ubd(X)) / 3
    + -(2 * Nf ** 6 * db(X) * L(X) ** 3 * Ld(X) * ubd(X)) / 3
    + 2 * Nf ** 6 * db(X) ** 2 * ebd(X) * L(X) * Q(X) * ubd(X)
    + -5 * Nf ** 4 * D * ebd(X) * H(X) ** 2 * L(X) * Qd(X) * ubd(X)
    + 7 * Nf ** 4 * D ** 2 * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + -Nf ** 4 * B(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + Nf ** 4 * Bb(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + Nf ** 4 * G(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + -Nf ** 4 * Gb(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + (Nf ** 3 * H(X) ** 2 * Hd(X) * L(X) ** 2 * Qd(X) * ubd(X)) / 2
    + -(3 * Nf ** 4 * H(X) ** 2 * Hd(X) * L(X) ** 2 * Qd(X) * ubd(X)) / 2
    + (Nf ** 4 * eb(X) * L(X) ** 3 * Qd(X) * ubd(X)) / 3
    + -(2 * Nf ** 6 * eb(X) * L(X) ** 3 * Qd(X) * ubd(X)) / 3
    + 4 * Nf ** 6 * db(X) * L(X) ** 2 * Q(X) * Qd(X) * ubd(X)
    - -(Nf ** 3 * D * H(X) ** 2 * L(X) ** 2 * ub(X) * ubd(X)) / 2
    + (5 * Nf ** 4 * D * H(X) ** 2 * L(X) ** 2 * ub(X) * ubd(X)) / 2
    + -(Nf ** 3 * db(X) ** 2 * ebd(X) ** 2 * ubd(X) ** 2) / 2
    + (Nf ** 6 * db(X) ** 2 * ebd(X) ** 2 * ubd(X) ** 2) / 2
    + -2 * Nf ** 6 * db(X) * ebd(X) * L(X) * Qd(X) * ubd(X) ** 2
    + Nf ** 3 * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2
    + Nf ** 6 * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2
    + -Nf ** 6 * db(X) * L(X) ** 2 * ub(X) * ubd(X) ** 2
    + 2 * Nf ** 2 * D * ebd(X) * H(X) ** 3 * L(X) * W(X)
    + 5 * Nf ** 2 * D ** 2 * H(X) ** 2 * L(X) ** 2 * W(X)
    + -2 * Nf ** 2 * B(X) * H(X) ** 2 * L(X) ** 2 * W(X)
    - (Nf * H(X) ** 3 * Hd(X) * L(X) ** 2 * W(X)) / 2
    + (3 * Nf ** 2 * H(X) ** 3 * Hd(X) * L(X) ** 2 * W(X)) / 2
    - -(Nf ** 3 * eb(X) * H(X) * L(X) ** 3 * W(X)) / 2
    + (3 * Nf ** 4 * eb(X) * H(X) * L(X) ** 3 * W(X)) / 2
    - (Nf ** 3 * db(X) * H(X) * L(X) ** 2 * Q(X) * W(X)) / 2
    + -(9 * Nf ** 4 * db(X) * H(X) * L(X) ** 2 * Q(X) * W(X)) / 2
    + Nf ** 4 * db(X) * ebd(X) * H(X) * L(X) * ubd(X) * W(X)
    + -2 * Nf ** 4 * D * db(X) * L(X) ** 2 * ubd(X) * W(X)
    - (Nf ** 3 * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * W(X)) / 2
    + -(3 * Nf ** 4 * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * W(X)) / 2
    + Nf * H(X) ** 2 * L(X) ** 2 * W(X) ** 2
    + 2 * Nf ** 2 * H(X) ** 2 * L(X) ** 2 * W(X) ** 2
    + -2 * Nf ** 2 * D * ebd(X) * H(X) ** 3 * L(X) * Wb(X)
    + 3 * Nf ** 2 * D ** 2 * H(X) ** 2 * L(X) ** 2 * Wb(X)
    + Nf ** 2 * Bb(X) * H(X) ** 2 * L(X) ** 2 * Wb(X)
    + -Nf ** 4 * db(X) * ebd(X) * H(X) * L(X) * ubd(X) * Wb(X)
    - (Nf ** 3 * D * db(X) * L(X) ** 2 * ubd(X) * Wb(X)) / 2
    + -(3 * Nf ** 4 * D * db(X) * L(X) ** 2 * ubd(X) * Wb(X)) / 2
    + (Nf ** 3 * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * Wb(X)) / 2
    + -(3 * Nf ** 4 * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * Wb(X)) / 2
    + Nf * H(X) ** 2 * L(X) ** 2 * Wb(X) ** 2
    + Nf ** 2 * H(X) ** 2 * L(X) ** 2 * Wb(X) ** 2
)

# No derivatives here
H11_LNV = (
    (Nf ** 3 * db(X) ** 2 * dbd(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    + (Nf ** 6 * db(X) ** 2 * dbd(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    + -Nf ** 6 * db(X) * dbd(X) * eb(X) * ebd(X) * H(X) ** 2 * L(X) ** 2
    + (Nf ** 3 * eb(X) ** 2 * ebd(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 4
    + -(Nf ** 4 * eb(X) ** 2 * ebd(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 4
    + (Nf ** 5 * eb(X) ** 2 * ebd(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 4
    + -(Nf ** 6 * eb(X) ** 2 * ebd(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 4
    + (Nf * H(X) ** 5 * Hd(X) ** 3 * L(X) ** 2) / 2
    + -(Nf ** 2 * H(X) ** 5 * Hd(X) ** 3 * L(X) ** 2) / 2
    + Nf ** 4 * eb(X) * H(X) ** 3 * Hd(X) ** 2 * L(X) ** 3
    + (Nf ** 2 * eb(X) ** 2 * H(X) * Hd(X) * L(X) ** 4) / 8
    + -(11 * Nf ** 3 * eb(X) ** 2 * H(X) * Hd(X) * L(X) ** 4) / 48
    + (17 * Nf ** 4 * eb(X) ** 2 * H(X) * Hd(X) * L(X) ** 4) / 48
    - -(11 * Nf ** 5 * eb(X) ** 2 * H(X) * Hd(X) * L(X) ** 4) / 48
    + (25 * Nf ** 6 * eb(X) ** 2 * H(X) * Hd(X) * L(X) ** 4) / 48
    + -(Nf ** 3 * ebd(X) * H(X) ** 4 * Hd(X) * L(X) ** 2 * Ld(X)) / 2
    + (Nf ** 4 * ebd(X) * H(X) ** 4 * Hd(X) * L(X) ** 2 * Ld(X)) / 2
    + -Nf ** 6 * db(X) * dbd(X) * H(X) ** 2 * L(X) ** 3 * Ld(X)
    + Nf ** 6 * eb(X) * ebd(X) * H(X) ** 2 * L(X) ** 3 * Ld(X)
    + -(3 * Nf ** 3 * H(X) ** 2 * L(X) ** 4 * Ld(X) ** 2) / 8
    + (Nf ** 4 * H(X) ** 2 * L(X) ** 4 * Ld(X) ** 2) / 8
    + (Nf ** 5 * H(X) ** 2 * L(X) ** 4 * Ld(X) ** 2) / 8
    + -(3 * Nf ** 6 * H(X) ** 2 * L(X) ** 4 * Ld(X) ** 2) / 8
    + 2 * Nf ** 6 * db(X) ** 2 * dbd(X) * ebd(X) * H(X) ** 2 * L(X) * Q(X)
    + -Nf ** 5 * db(X) * eb(X) * ebd(X) ** 2 * H(X) ** 2 * L(X) * Q(X)
    + Nf ** 6 * db(X) * eb(X) * ebd(X) ** 2 * H(X) ** 2 * L(X) * Q(X)
    + -3 * Nf ** 4 * db(X) * H(X) ** 3 * Hd(X) ** 2 * L(X) ** 2 * Q(X)
    + (Nf ** 4 * db(X) * eb(X) * H(X) * Hd(X) * L(X) ** 3 * Q(X)) / 3
    - -(Nf ** 5 * db(X) * eb(X) * H(X) * Hd(X) * L(X) ** 3 * Q(X)) / 2
    + (25 * Nf ** 6 * db(X) * eb(X) * H(X) * Hd(X) * L(X) ** 3 * Q(X)) / 6
    + -3 * Nf ** 6 * db(X) * ebd(X) * H(X) ** 2 * L(X) ** 2 * Ld(X) * Q(X)
    + (Nf ** 3 * db(X) ** 2 * ebd(X) ** 2 * H(X) ** 2 * Q(X) ** 2) / 2
    + -(Nf ** 4 * db(X) ** 2 * ebd(X) ** 2 * H(X) ** 2 * Q(X) ** 2) / 2
    + (Nf ** 5 * db(X) ** 2 * ebd(X) ** 2 * H(X) ** 2 * Q(X) ** 2) / 2
    + -(Nf ** 6 * db(X) ** 2 * ebd(X) ** 2 * H(X) ** 2 * Q(X) ** 2) / 2
    + (3 * Nf ** 3 * db(X) ** 2 * H(X) * Hd(X) * L(X) ** 2 * Q(X) ** 2) / 4
    + -(Nf ** 4 * db(X) ** 2 * H(X) * Hd(X) * L(X) ** 2 * Q(X) ** 2) / 4
    - (Nf ** 5 * db(X) ** 2 * H(X) * Hd(X) * L(X) ** 2 * Q(X) ** 2) / 4
    + -(25 * Nf ** 6 * db(X) ** 2 * H(X) * Hd(X) * L(X) ** 2 * Q(X) ** 2) / 4
    + (Nf ** 3 * dbd(X) * H(X) ** 4 * Hd(X) * L(X) ** 2 * Qd(X)) / 2
    + -(Nf ** 4 * dbd(X) * H(X) ** 4 * Hd(X) * L(X) ** 2 * Qd(X)) / 2
    + Nf ** 6 * dbd(X) * eb(X) * H(X) ** 2 * L(X) ** 3 * Qd(X)
    + -Nf ** 4 * ebd(X) * H(X) ** 4 * Hd(X) * L(X) * Q(X) * Qd(X)
    + 6 * Nf ** 6 * db(X) * dbd(X) * H(X) ** 2 * L(X) ** 2 * Q(X) * Qd(X)
    + -3 * Nf ** 6 * eb(X) * ebd(X) * H(X) ** 2 * L(X) ** 2 * Q(X) * Qd(X)
    + 3 * Nf ** 6 * H(X) ** 2 * L(X) ** 3 * Ld(X) * Q(X) * Qd(X)
    + -6 * Nf ** 6 * db(X) * ebd(X) * H(X) ** 2 * L(X) * Q(X) ** 2 * Qd(X)
    + (3 * Nf ** 3 * H(X) ** 2 * L(X) ** 2 * Q(X) ** 2 * Qd(X) ** 2) / 2
    + -(9 * Nf ** 6 * H(X) ** 2 * L(X) ** 2 * Q(X) ** 2 * Qd(X) ** 2) / 2
    + Nf ** 4 * H(X) ** 4 * Hd(X) * L(X) ** 2 * Q(X) * ub(X)
    - -(Nf ** 5 * eb(X) * H(X) ** 2 * L(X) ** 3 * Q(X) * ub(X)) / 2
    + (5 * Nf ** 6 * eb(X) * H(X) ** 2 * L(X) ** 3 * Q(X) * ub(X)) / 2
    - -(Nf ** 5 * db(X) * H(X) ** 2 * L(X) ** 2 * Q(X) ** 2 * ub(X)) / 2
    + (15 * Nf ** 6 * db(X) * H(X) ** 2 * L(X) ** 2 * Q(X) ** 2 * ub(X)) / 2
    + -Nf ** 4 * db(X) * ebd(X) * H(X) ** 3 * Hd(X) ** 2 * L(X) * ubd(X)
    + 2 * Nf ** 6 * db(X) ** 2 * dbd(X) * H(X) * Hd(X) * L(X) ** 2 * ubd(X)
    + -2 * Nf ** 6 * db(X) * eb(X) * ebd(X) * H(X) * Hd(X) * L(X) ** 2 * ubd(X)
    + Nf ** 6 * db(X) * ebd(X) ** 2 * H(X) ** 2 * L(X) * Ld(X) * ubd(X)
    + -(Nf ** 4 * db(X) * H(X) * Hd(X) * L(X) ** 3 * Ld(X) * ubd(X)) / 3
    + (5 * Nf ** 6 * db(X) * H(X) * Hd(X) * L(X) ** 3 * Ld(X) * ubd(X)) / 3
    + -4 * Nf ** 6 * db(X) ** 2 * ebd(X) * H(X) * Hd(X) * L(X) * Q(X) * ubd(X)
    + 4 * Nf ** 6 * db(X) * dbd(X) * ebd(X) * H(X) ** 2 * L(X) * Qd(X) * ubd(X)
    + -Nf ** 6 * eb(X) * ebd(X) ** 2 * H(X) ** 2 * L(X) * Qd(X) * ubd(X)
    + (Nf ** 3 * H(X) ** 3 * Hd(X) ** 2 * L(X) ** 2 * Qd(X) * ubd(X)) / 2
    + -(3 * Nf ** 4 * H(X) ** 3 * Hd(X) ** 2 * L(X) ** 2 * Qd(X) * ubd(X)) / 2
    + (Nf ** 4 * eb(X) * H(X) * Hd(X) * L(X) ** 3 * Qd(X) * ubd(X)) / 3
    + -(5 * Nf ** 6 * eb(X) * H(X) * Hd(X) * L(X) ** 3 * Qd(X) * ubd(X)) / 3
    + Nf ** 5 * ebd(X) * H(X) ** 2 * L(X) ** 2 * Ld(X) * Qd(X) * ubd(X)
    + -3 * Nf ** 6 * ebd(X) * H(X) ** 2 * L(X) ** 2 * Ld(X) * Qd(X) * ubd(X)
    + 2 * Nf ** 6 * db(X) * ebd(X) ** 2 * H(X) ** 2 * Q(X) * Qd(X) * ubd(X)
    + -10 * Nf ** 6 * db(X) * H(X) * Hd(X) * L(X) ** 2 * Q(X) * Qd(X) * ubd(X)
    + Nf ** 5 * dbd(X) * H(X) ** 2 * L(X) ** 2 * Qd(X) ** 2 * ubd(X)
    + -3 * Nf ** 6 * dbd(X) * H(X) ** 2 * L(X) ** 2 * Qd(X) ** 2 * ubd(X)
    + 6 * Nf ** 6 * ebd(X) * H(X) ** 2 * L(X) * Q(X) * Qd(X) ** 2 * ubd(X)
    + -2 * Nf ** 6 * db(X) * dbd(X) * H(X) ** 2 * L(X) ** 2 * ub(X) * ubd(X)
    + Nf ** 6 * eb(X) * ebd(X) * H(X) ** 2 * L(X) ** 2 * ub(X) * ubd(X)
    + -Nf ** 6 * H(X) ** 2 * L(X) ** 3 * Ld(X) * ub(X) * ubd(X)
    + 4 * Nf ** 6 * db(X) * ebd(X) * H(X) ** 2 * L(X) * Q(X) * ub(X) * ubd(X)
    + -6 * Nf ** 6 * H(X) ** 2 * L(X) ** 2 * Q(X) * Qd(X) * ub(X) * ubd(X)
    + (Nf ** 3 * db(X) ** 2 * ebd(X) ** 2 * H(X) * Hd(X) * ubd(X) ** 2) / 2
    + -(Nf ** 6 * db(X) ** 2 * ebd(X) ** 2 * H(X) * Hd(X) * ubd(X) ** 2) / 2
    + (Nf ** 3 * db(X) ** 2 * Hd(X) ** 2 * L(X) ** 2 * ubd(X) ** 2) / 2
    + -(Nf ** 6 * db(X) ** 2 * Hd(X) ** 2 * L(X) ** 2 * ubd(X) ** 2) / 2
    + 4 * Nf ** 6 * db(X) * ebd(X) * H(X) * Hd(X) * L(X) * Qd(X) * ubd(X) ** 2
    + -(3 * Nf ** 3 * ebd(X) ** 2 * H(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2) / 4
    + (Nf ** 4 * ebd(X) ** 2 * H(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2) / 4
    - -(Nf ** 5 * ebd(X) ** 2 * H(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2) / 4
    + (5 * Nf ** 6 * ebd(X) ** 2 * H(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2) / 4
    + -(Nf ** 3 * H(X) * Hd(X) * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2) / 2
    + (Nf ** 4 * H(X) * Hd(X) * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2) / 2
    + -(Nf ** 5 * H(X) * Hd(X) * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2) / 2
    + (5 * Nf ** 6 * H(X) * Hd(X) * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2) / 2
    + -2 * Nf ** 6 * db(X) * H(X) * Hd(X) * L(X) ** 2 * ub(X) * ubd(X) ** 2
    + 2 * Nf ** 6 * ebd(X) * H(X) ** 2 * L(X) * Qd(X) * ub(X) * ubd(X) ** 2
    + -(Nf ** 3 * H(X) ** 2 * L(X) ** 2 * ub(X) ** 2 * ubd(X) ** 2) / 2
    + (Nf ** 6 * H(X) ** 2 * L(X) ** 2 * ub(X) ** 2 * ubd(X) ** 2) / 2
)
