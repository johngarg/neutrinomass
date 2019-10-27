#!/usr/bin/env python3

"""Results from the Hilbert Series for the SMEFT up to dimension 7 and for the
ΔL = 2 SMEFT up to dimension 11."""

from sympy.abc import X
from sympy import symbols, Function

Nf, D = symbols("Nf D")
H, Hd, L, Ld, Q, Qd = symbols("H Hd L Ld Q Qd", cls=Function)
eb, ebd, ub, ubd, db, dbd = symbols("eb ebd ub ubd db dbd", cls=Function)
G, Gb, W, Wb, B, Bb = symbols("G Gb W Wb B Bb", cls=Function)

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