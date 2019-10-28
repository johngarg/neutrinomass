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

H11_LNV = (
    (3 * Nf * D ** 4 * ebd(X) ** 2 * H(X) ** 4) / 2
    + (7 * Nf ** 2 * D ** 4 * ebd(X) ** 2 * H(X) ** 4) / 2
    - (Nf * D ** 2 * B(X) * ebd(X) ** 2 * H(X) ** 4) / 2
    + -(Nf ** 2 * D ** 2 * B(X) * ebd(X) ** 2 * H(X) ** 4) / 2
    - (Nf * D ** 2 * Bb(X) * ebd(X) ** 2 * H(X) ** 4) / 2
    + -(Nf ** 2 * D ** 2 * Bb(X) * ebd(X) ** 2 * H(X) ** 4) / 2
    + (Nf * D ** 2 * ebd(X) ** 2 * H(X) ** 5 * Hd(X)) / 2
    + -(Nf ** 2 * D ** 2 * ebd(X) ** 2 * H(X) ** 5 * Hd(X)) / 2
    + 10 * Nf ** 2 * D ** 5 * ebd(X) * H(X) ** 3 * L(X)
    + 11 * Nf ** 2 * D ** 3 * B(X) * ebd(X) * H(X) ** 3 * L(X)
    + -Nf ** 2 * D * B(X) ** 2 * ebd(X) * H(X) ** 3 * L(X)
    + 11 * Nf ** 2 * D ** 3 * Bb(X) * ebd(X) * H(X) ** 3 * L(X)
    + -Nf ** 2 * D * B(X) * Bb(X) * ebd(X) * H(X) ** 3 * L(X)
    + Nf ** 2 * D * Bb(X) ** 2 * ebd(X) * H(X) ** 3 * L(X)
    + -9 * Nf ** 4 * D ** 2 * db(X) * dbd(X) * ebd(X) * H(X) ** 3 * L(X)
    - (Nf ** 3 * D ** 2 * eb(X) * ebd(X) ** 2 * H(X) ** 3 * L(X)) / 2
    + -(9 * Nf ** 4 * D ** 2 * eb(X) * ebd(X) ** 2 * H(X) ** 3 * L(X)) / 2
    + Nf ** 2 * D * ebd(X) * G(X) ** 2 * H(X) ** 3 * L(X)
    + -Nf ** 2 * D * ebd(X) * G(X) * Gb(X) * H(X) ** 3 * L(X)
    + Nf ** 2 * D * ebd(X) * Gb(X) ** 2 * H(X) ** 3 * L(X)
    + -14 * Nf ** 2 * D ** 3 * ebd(X) * H(X) ** 4 * Hd(X) * L(X)
    + Nf ** 2 * D * B(X) * ebd(X) * H(X) ** 4 * Hd(X) * L(X)
    + -Nf ** 2 * D * Bb(X) * ebd(X) * H(X) ** 4 * Hd(X) * L(X)
    + Nf ** 2 * D * ebd(X) * H(X) ** 5 * Hd(X) ** 2 * L(X)
    + 2 * Nf * D ** 6 * H(X) ** 2 * L(X) ** 2
    + -2 * Nf ** 2 * D ** 6 * H(X) ** 2 * L(X) ** 2
    - 2 * Nf * D ** 4 * B(X) * H(X) ** 2 * L(X) ** 2
    + 12 * Nf ** 2 * D ** 4 * B(X) * H(X) ** 2 * L(X) ** 2
    + -(5 * Nf * D ** 2 * B(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    + (13 * Nf ** 2 * D ** 2 * B(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    - -(Nf * B(X) ** 3 * H(X) ** 2 * L(X) ** 2) / 2
    + (Nf ** 2 * B(X) ** 3 * H(X) ** 2 * L(X) ** 2) / 2
    - Nf * D ** 4 * Bb(X) * H(X) ** 2 * L(X) ** 2
    + -9 * Nf ** 2 * D ** 4 * Bb(X) * H(X) ** 2 * L(X) ** 2
    + Nf * D ** 2 * B(X) * Bb(X) * H(X) ** 2 * L(X) ** 2
    + 6 * Nf ** 2 * D ** 2 * B(X) * Bb(X) * H(X) ** 2 * L(X) ** 2
    + -(3 * Nf * D ** 2 * Bb(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    + (7 * Nf ** 2 * D ** 2 * Bb(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    - -(Nf * B(X) * Bb(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    + (Nf ** 2 * B(X) * Bb(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    - -(Nf ** 3 * D ** 3 * db(X) * dbd(X) * H(X) ** 2 * L(X) ** 2) / 2
    + (45 * Nf ** 4 * D ** 3 * db(X) * dbd(X) * H(X) ** 2 * L(X) ** 2) / 2
    + -5 * Nf ** 4 * D * B(X) * db(X) * dbd(X) * H(X) ** 2 * L(X) ** 2
    - (Nf ** 3 * D * Bb(X) * db(X) * dbd(X) * H(X) ** 2 * L(X) ** 2) / 2
    + -(7 * Nf ** 4 * D * Bb(X) * db(X) * dbd(X) * H(X) ** 2 * L(X) ** 2) / 2
    + (Nf ** 3 * db(X) ** 2 * dbd(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    + -(Nf ** 6 * db(X) ** 2 * dbd(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    - (Nf ** 3 * D ** 3 * eb(X) * ebd(X) * H(X) ** 2 * L(X) ** 2) / 2
    + -(45 * Nf ** 4 * D ** 3 * eb(X) * ebd(X) * H(X) ** 2 * L(X) ** 2) / 2
    + 5 * Nf ** 4 * D * B(X) * eb(X) * ebd(X) * H(X) ** 2 * L(X) ** 2
    - -(Nf ** 3 * D * Bb(X) * eb(X) * ebd(X) * H(X) ** 2 * L(X) ** 2) / 2
    + (7 * Nf ** 4 * D * Bb(X) * eb(X) * ebd(X) * H(X) ** 2 * L(X) ** 2) / 2
    + -Nf ** 6 * db(X) * dbd(X) * eb(X) * ebd(X) * H(X) ** 2 * L(X) ** 2
    + (Nf ** 3 * eb(X) ** 2 * ebd(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 4
    + -(Nf ** 4 * eb(X) ** 2 * ebd(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 4
    + (Nf ** 5 * eb(X) ** 2 * ebd(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 4
    + -(Nf ** 6 * eb(X) ** 2 * ebd(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 4
    + 5 * Nf ** 4 * D * db(X) * dbd(X) * G(X) * H(X) ** 2 * L(X) ** 2
    + -(5 * Nf * D ** 2 * G(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    + (13 * Nf ** 2 * D ** 2 * G(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    - -Nf * B(X) * G(X) ** 2 * H(X) ** 2 * L(X) ** 2
    + Nf ** 2 * B(X) * G(X) ** 2 * H(X) ** 2 * L(X) ** 2
    + Nf ** 2 * G(X) ** 3 * H(X) ** 2 * L(X) ** 2
    - -(Nf ** 3 * D * db(X) * dbd(X) * Gb(X) * H(X) ** 2 * L(X) ** 2) / 2
    + (7 * Nf ** 4 * D * db(X) * dbd(X) * Gb(X) * H(X) ** 2 * L(X) ** 2) / 2
    + -Nf * D ** 2 * G(X) * Gb(X) * H(X) ** 2 * L(X) ** 2
    + 6 * Nf ** 2 * D ** 2 * G(X) * Gb(X) * H(X) ** 2 * L(X) ** 2
    - -(Nf * Bb(X) * G(X) * Gb(X) * H(X) ** 2 * L(X) ** 2) / 2
    + (Nf ** 2 * Bb(X) * G(X) * Gb(X) * H(X) ** 2 * L(X) ** 2) / 2
    + -(3 * Nf * D ** 2 * Gb(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    + (7 * Nf ** 2 * D ** 2 * Gb(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    - -(Nf * B(X) * Gb(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    + (Nf ** 2 * B(X) * Gb(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    - -(Nf * G(X) * Gb(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    + (Nf ** 2 * G(X) * Gb(X) ** 2 * H(X) ** 2 * L(X) ** 2) / 2
    + (Nf * Gb(X) ** 3 * H(X) ** 2 * L(X) ** 2) / 2
    + -(Nf ** 2 * Gb(X) ** 3 * H(X) ** 2 * L(X) ** 2) / 2
    + (3 * Nf * D ** 4 * H(X) ** 3 * Hd(X) * L(X) ** 2) / 2
    + -(63 * Nf ** 2 * D ** 4 * H(X) ** 3 * Hd(X) * L(X) ** 2) / 2
    - Nf * D ** 2 * B(X) * H(X) ** 3 * Hd(X) * L(X) ** 2
    + -12 * Nf ** 2 * D ** 2 * B(X) * H(X) ** 3 * Hd(X) * L(X) ** 2
    + (Nf * B(X) ** 2 * H(X) ** 3 * Hd(X) * L(X) ** 2) / 2
    + -(Nf ** 2 * B(X) ** 2 * H(X) ** 3 * Hd(X) * L(X) ** 2) / 2
    - (3 * Nf * D ** 2 * Bb(X) * H(X) ** 3 * Hd(X) * L(X) ** 2) / 2
    + -(15 * Nf ** 2 * D ** 2 * Bb(X) * H(X) ** 3 * Hd(X) * L(X) ** 2) / 2
    + (Nf * Bb(X) ** 2 * H(X) ** 3 * Hd(X) * L(X) ** 2) / 2
    + -(Nf ** 2 * Bb(X) ** 2 * H(X) ** 3 * Hd(X) * L(X) ** 2) / 2
    - (Nf ** 3 * D * db(X) * dbd(X) * H(X) ** 3 * Hd(X) * L(X) ** 2) / 2
    + -(9 * Nf ** 4 * D * db(X) * dbd(X) * H(X) ** 3 * Hd(X) * L(X) ** 2) / 2
    - (Nf ** 3 * D * eb(X) * ebd(X) * H(X) ** 3 * Hd(X) * L(X) ** 2) / 2
    + -(9 * Nf ** 4 * D * eb(X) * ebd(X) * H(X) ** 3 * Hd(X) * L(X) ** 2) / 2
    + (Nf * G(X) ** 2 * H(X) ** 3 * Hd(X) * L(X) ** 2) / 2
    + -(Nf ** 2 * G(X) ** 2 * H(X) ** 3 * Hd(X) * L(X) ** 2) / 2
    + (Nf * Gb(X) ** 2 * H(X) ** 3 * Hd(X) * L(X) ** 2) / 2
    + -(Nf ** 2 * Gb(X) ** 2 * H(X) ** 3 * Hd(X) * L(X) ** 2) / 2
    + (3 * Nf * D ** 2 * H(X) ** 4 * Hd(X) ** 2 * L(X) ** 2) / 2
    + -(15 * Nf ** 2 * D ** 2 * H(X) ** 4 * Hd(X) ** 2 * L(X) ** 2) / 2
    - (Nf * B(X) * H(X) ** 4 * Hd(X) ** 2 * L(X) ** 2) / 2
    + -(Nf ** 2 * B(X) * H(X) ** 4 * Hd(X) ** 2 * L(X) ** 2) / 2
    + (Nf * H(X) ** 5 * Hd(X) ** 3 * L(X) ** 2) / 2
    + (Nf ** 2 * H(X) ** 5 * Hd(X) ** 3 * L(X) ** 2) / 2
    + -10 * Nf ** 4 * D ** 4 * eb(X) * H(X) * L(X) ** 3
    + 11 * Nf ** 4 * D ** 2 * B(X) * eb(X) * H(X) * L(X) ** 3
    + Nf ** 4 * B(X) ** 2 * eb(X) * H(X) * L(X) ** 3
    - -(Nf ** 2 * D ** 2 * Bb(X) * eb(X) * H(X) * L(X) ** 3) / 3
    + (16 * Nf ** 4 * D ** 2 * Bb(X) * eb(X) * H(X) * L(X) ** 3) / 3
    + -(Nf ** 2 * Bb(X) ** 2 * eb(X) * H(X) * L(X) ** 3) / 3
    + (2 * Nf ** 4 * Bb(X) ** 2 * eb(X) * H(X) * L(X) ** 3) / 3
    + -5 * Nf ** 6 * D * db(X) * dbd(X) * eb(X) * H(X) * L(X) ** 3
    - (Nf ** 5 * D * eb(X) ** 2 * ebd(X) * H(X) * L(X) ** 3) / 2
    + -(5 * Nf ** 6 * D * eb(X) ** 2 * ebd(X) * H(X) * L(X) ** 3) / 2
    + Nf ** 4 * eb(X) * G(X) ** 2 * H(X) * L(X) ** 3
    + -(Nf ** 2 * eb(X) * Gb(X) ** 2 * H(X) * L(X) ** 3) / 3
    + (2 * Nf ** 4 * eb(X) * Gb(X) ** 2 * H(X) * L(X) ** 3) / 3
    - -Nf ** 3 * D ** 2 * eb(X) * H(X) ** 2 * Hd(X) * L(X) ** 3
    + 17 * Nf ** 4 * D ** 2 * eb(X) * H(X) ** 2 * Hd(X) * L(X) ** 3
    - -(Nf ** 3 * B(X) * eb(X) * H(X) ** 2 * Hd(X) * L(X) ** 3) / 2
    + (3 * Nf ** 4 * B(X) * eb(X) * H(X) ** 2 * Hd(X) * L(X) ** 3) / 2
    + -Nf ** 4 * eb(X) * H(X) ** 3 * Hd(X) ** 2 * L(X) ** 3
    + (9 * Nf ** 3 * D ** 2 * eb(X) ** 2 * L(X) ** 4) / 8
    + (Nf ** 4 * D ** 2 * eb(X) ** 2 * L(X) ** 4) / 8
    - -(Nf ** 5 * D ** 2 * eb(X) ** 2 * L(X) ** 4) / 8
    + (15 * Nf ** 6 * D ** 2 * eb(X) ** 2 * L(X) ** 4) / 8
    - (3 * Nf ** 3 * B(X) * eb(X) ** 2 * L(X) ** 4) / 8
    + -(Nf ** 4 * B(X) * eb(X) ** 2 * L(X) ** 4) / 8
    - (Nf ** 5 * B(X) * eb(X) ** 2 * L(X) ** 4) / 8
    + (3 * Nf ** 6 * B(X) * eb(X) ** 2 * L(X) ** 4) / 8
    + -(Nf ** 2 * eb(X) ** 2 * H(X) * Hd(X) * L(X) ** 4) / 8
    + (11 * Nf ** 3 * eb(X) ** 2 * H(X) * Hd(X) * L(X) ** 4) / 48
    + -(17 * Nf ** 4 * eb(X) ** 2 * H(X) * Hd(X) * L(X) ** 4) / 48
    - (11 * Nf ** 5 * eb(X) ** 2 * H(X) * Hd(X) * L(X) ** 4) / 48
    + -(25 * Nf ** 6 * eb(X) ** 2 * H(X) * Hd(X) * L(X) ** 4) / 48
    + Nf ** 4 * D * ebd(X) ** 2 * H(X) ** 4 * L(X) * Ld(X)
    - -(Nf ** 3 * D ** 2 * ebd(X) * H(X) ** 3 * L(X) ** 2 * Ld(X)) / 2
    + (27 * Nf ** 4 * D ** 2 * ebd(X) * H(X) ** 3 * L(X) ** 2 * Ld(X)) / 2
    - -(Nf ** 3 * B(X) * ebd(X) * H(X) ** 3 * L(X) ** 2 * Ld(X)) / 2
    + (Nf ** 4 * B(X) * ebd(X) * H(X) ** 3 * L(X) ** 2 * Ld(X)) / 2
    + -(Nf ** 3 * Bb(X) * ebd(X) * H(X) ** 3 * L(X) ** 2 * Ld(X)) / 2
    + (Nf ** 4 * Bb(X) * ebd(X) * H(X) ** 3 * L(X) ** 2 * Ld(X)) / 2
    + -(Nf ** 3 * ebd(X) * H(X) ** 4 * Hd(X) * L(X) ** 2 * Ld(X)) / 2
    + (Nf ** 4 * ebd(X) * H(X) ** 4 * Hd(X) * L(X) ** 2 * Ld(X)) / 2
    - -Nf ** 3 * D ** 3 * H(X) ** 2 * L(X) ** 3 * Ld(X)
    + 19 * Nf ** 4 * D ** 3 * H(X) ** 2 * L(X) ** 3 * Ld(X)
    - (Nf ** 3 * D * B(X) * H(X) ** 2 * L(X) ** 3 * Ld(X)) / 2
    + -(9 * Nf ** 4 * D * B(X) * H(X) ** 2 * L(X) ** 3 * Ld(X)) / 2
    + (Nf ** 2 * D * Bb(X) * H(X) ** 2 * L(X) ** 3 * Ld(X)) / 3
    - -(Nf ** 3 * D * Bb(X) * H(X) ** 2 * L(X) ** 3 * Ld(X)) / 2
    + (19 * Nf ** 4 * D * Bb(X) * H(X) ** 2 * L(X) ** 3 * Ld(X)) / 6
    + -Nf ** 6 * db(X) * dbd(X) * H(X) ** 2 * L(X) ** 3 * Ld(X)
    + Nf ** 6 * eb(X) * ebd(X) * H(X) ** 2 * L(X) ** 3 * Ld(X)
    - -Nf ** 3 * D * H(X) ** 3 * Hd(X) * L(X) ** 3 * Ld(X)
    + 5 * Nf ** 4 * D * H(X) ** 3 * Hd(X) * L(X) ** 3 * Ld(X)
    - (Nf ** 3 * D * eb(X) * H(X) * L(X) ** 4 * Ld(X)) / 4
    - -(Nf ** 4 * D * eb(X) * H(X) * L(X) ** 4 * Ld(X)) / 8
    - (3 * Nf ** 5 * D * eb(X) * H(X) * L(X) ** 4 * Ld(X)) / 4
    + -(25 * Nf ** 6 * D * eb(X) * H(X) * L(X) ** 4 * Ld(X)) / 8
    + (3 * Nf ** 3 * H(X) ** 2 * L(X) ** 4 * Ld(X) ** 2) / 8
    + -(Nf ** 4 * H(X) ** 2 * L(X) ** 4 * Ld(X) ** 2) / 8
    + (Nf ** 5 * H(X) ** 2 * L(X) ** 4 * Ld(X) ** 2) / 8
    + (3 * Nf ** 6 * H(X) ** 2 * L(X) ** 4 * Ld(X) ** 2) / 8
    - -(Nf ** 3 * D ** 2 * db(X) * ebd(X) ** 2 * H(X) ** 3 * Q(X)) / 2
    + (9 * Nf ** 4 * D ** 2 * db(X) * ebd(X) ** 2 * H(X) ** 3 * Q(X)) / 2
    + -45 * Nf ** 4 * D ** 3 * db(X) * ebd(X) * H(X) ** 2 * L(X) * Q(X)
    + 10 * Nf ** 4 * D * B(X) * db(X) * ebd(X) * H(X) ** 2 * L(X) * Q(X)
    + -7 * Nf ** 4 * D * Bb(X) * db(X) * ebd(X) * H(X) ** 2 * L(X) * Q(X)
    + 2 * Nf ** 6 * db(X) ** 2 * dbd(X) * ebd(X) * H(X) ** 2 * L(X) * Q(X)
    + -Nf ** 5 * db(X) * eb(X) * ebd(X) ** 2 * H(X) ** 2 * L(X) * Q(X)
    + Nf ** 6 * db(X) * eb(X) * ebd(X) ** 2 * H(X) ** 2 * L(X) * Q(X)
    + -10 * Nf ** 4 * D * db(X) * ebd(X) * G(X) * H(X) ** 2 * L(X) * Q(X)
    + 7 * Nf ** 4 * D * db(X) * ebd(X) * Gb(X) * H(X) ** 2 * L(X) * Q(X)
    + -9 * Nf ** 4 * D * db(X) * ebd(X) * H(X) ** 3 * Hd(X) * L(X) * Q(X)
    + 30 * Nf ** 4 * D ** 4 * db(X) * H(X) * L(X) ** 2 * Q(X)
    + -33 * Nf ** 4 * D ** 2 * B(X) * db(X) * H(X) * L(X) ** 2 * Q(X)
    + 3 * Nf ** 4 * B(X) ** 2 * db(X) * H(X) * L(X) ** 2 * Q(X)
    + -16 * Nf ** 4 * D ** 2 * Bb(X) * db(X) * H(X) * L(X) ** 2 * Q(X)
    + 2 * Nf ** 4 * Bb(X) ** 2 * db(X) * H(X) * L(X) ** 2 * Q(X)
    + -15 * Nf ** 6 * D * db(X) ** 2 * dbd(X) * H(X) * L(X) ** 2 * Q(X)
    + 15 * Nf ** 6 * D * db(X) * eb(X) * ebd(X) * H(X) * L(X) ** 2 * Q(X)
    + -33 * Nf ** 4 * D ** 2 * db(X) * G(X) * H(X) * L(X) ** 2 * Q(X)
    + 6 * Nf ** 4 * B(X) * db(X) * G(X) * H(X) * L(X) ** 2 * Q(X)
    + -9 * Nf ** 4 * db(X) * G(X) ** 2 * H(X) * L(X) ** 2 * Q(X)
    + 16 * Nf ** 4 * D ** 2 * db(X) * Gb(X) * H(X) * L(X) ** 2 * Q(X)
    + -2 * Nf ** 4 * Bb(X) * db(X) * Gb(X) * H(X) * L(X) ** 2 * Q(X)
    + 4 * Nf ** 4 * db(X) * Gb(X) ** 2 * H(X) * L(X) ** 2 * Q(X)
    - -Nf ** 3 * D ** 2 * db(X) * H(X) ** 2 * Hd(X) * L(X) ** 2 * Q(X)
    + 51 * Nf ** 4 * D ** 2 * db(X) * H(X) ** 2 * Hd(X) * L(X) ** 2 * Q(X)
    - -(Nf ** 3 * B(X) * db(X) * H(X) ** 2 * Hd(X) * L(X) ** 2 * Q(X)) / 2
    + (9 * Nf ** 4 * B(X) * db(X) * H(X) ** 2 * Hd(X) * L(X) ** 2 * Q(X)) / 2
    - -(Nf ** 3 * db(X) * G(X) * H(X) ** 2 * Hd(X) * L(X) ** 2 * Q(X)) / 2
    + (9 * Nf ** 4 * db(X) * G(X) * H(X) ** 2 * Hd(X) * L(X) ** 2 * Q(X)) / 2
    + -3 * Nf ** 4 * db(X) * H(X) ** 3 * Hd(X) ** 2 * L(X) ** 2 * Q(X)
    + 15 * Nf ** 6 * D ** 2 * db(X) * eb(X) * L(X) ** 3 * Q(X)
    + -3 * Nf ** 6 * B(X) * db(X) * eb(X) * L(X) ** 3 * Q(X)
    + 3 * Nf ** 6 * db(X) * eb(X) * G(X) * L(X) ** 3 * Q(X)
    + -(Nf ** 4 * db(X) * eb(X) * H(X) * Hd(X) * L(X) ** 3 * Q(X)) / 3
    - (Nf ** 5 * db(X) * eb(X) * H(X) * Hd(X) * L(X) ** 3 * Q(X)) / 2
    + -(25 * Nf ** 6 * db(X) * eb(X) * H(X) * Hd(X) * L(X) ** 3 * Q(X)) / 6
    + 3 * Nf ** 6 * db(X) * ebd(X) * H(X) ** 2 * L(X) ** 2 * Ld(X) * Q(X)
    - -(3 * Nf ** 5 * D * db(X) * H(X) * L(X) ** 3 * Ld(X) * Q(X)) / 2
    + (25 * Nf ** 6 * D * db(X) * H(X) * L(X) ** 3 * Ld(X) * Q(X)) / 2
    + -(Nf ** 3 * db(X) ** 2 * ebd(X) ** 2 * H(X) ** 2 * Q(X) ** 2) / 2
    + (Nf ** 4 * db(X) ** 2 * ebd(X) ** 2 * H(X) ** 2 * Q(X) ** 2) / 2
    + -(Nf ** 5 * db(X) ** 2 * ebd(X) ** 2 * H(X) ** 2 * Q(X) ** 2) / 2
    + (Nf ** 6 * db(X) ** 2 * ebd(X) ** 2 * H(X) ** 2 * Q(X) ** 2) / 2
    + -15 * Nf ** 6 * D * db(X) ** 2 * ebd(X) * H(X) * L(X) * Q(X) ** 2
    + (9 * Nf ** 3 * D ** 2 * db(X) ** 2 * L(X) ** 2 * Q(X) ** 2) / 2
    + -(45 * Nf ** 6 * D ** 2 * db(X) ** 2 * L(X) ** 2 * Q(X) ** 2) / 2
    - (3 * Nf ** 3 * B(X) * db(X) ** 2 * L(X) ** 2 * Q(X) ** 2) / 2
    + -(9 * Nf ** 6 * B(X) * db(X) ** 2 * L(X) ** 2 * Q(X) ** 2) / 2
    + 9 * Nf ** 6 * db(X) ** 2 * G(X) * L(X) ** 2 * Q(X) ** 2
    + -(3 * Nf ** 3 * db(X) ** 2 * H(X) * Hd(X) * L(X) ** 2 * Q(X) ** 2) / 4
    + (Nf ** 4 * db(X) ** 2 * H(X) * Hd(X) * L(X) ** 2 * Q(X) ** 2) / 4
    - -(Nf ** 5 * db(X) ** 2 * H(X) * Hd(X) * L(X) ** 2 * Q(X) ** 2) / 4
    + (25 * Nf ** 6 * db(X) ** 2 * H(X) * Hd(X) * L(X) ** 2 * Q(X) ** 2) / 4
    + -2 * Nf ** 4 * D * dbd(X) * ebd(X) * H(X) ** 4 * L(X) * Qd(X)
    - (Nf ** 3 * D ** 2 * dbd(X) * H(X) ** 3 * L(X) ** 2 * Qd(X)) / 2
    + -(27 * Nf ** 4 * D ** 2 * dbd(X) * H(X) ** 3 * L(X) ** 2 * Qd(X)) / 2
    - (Nf ** 3 * B(X) * dbd(X) * H(X) ** 3 * L(X) ** 2 * Qd(X)) / 2
    + -(Nf ** 4 * B(X) * dbd(X) * H(X) ** 3 * L(X) ** 2 * Qd(X)) / 2
    + (Nf ** 3 * Bb(X) * dbd(X) * H(X) ** 3 * L(X) ** 2 * Qd(X)) / 2
    + -(Nf ** 4 * Bb(X) * dbd(X) * H(X) ** 3 * L(X) ** 2 * Qd(X)) / 2
    - (Nf ** 3 * dbd(X) * G(X) * H(X) ** 3 * L(X) ** 2 * Qd(X)) / 2
    + -(Nf ** 4 * dbd(X) * G(X) * H(X) ** 3 * L(X) ** 2 * Qd(X)) / 2
    + (Nf ** 3 * dbd(X) * Gb(X) * H(X) ** 3 * L(X) ** 2 * Qd(X)) / 2
    + -(Nf ** 4 * dbd(X) * Gb(X) * H(X) ** 3 * L(X) ** 2 * Qd(X)) / 2
    + (Nf ** 3 * dbd(X) * H(X) ** 4 * Hd(X) * L(X) ** 2 * Qd(X)) / 2
    + -(Nf ** 4 * dbd(X) * H(X) ** 4 * Hd(X) * L(X) ** 2 * Qd(X)) / 2
    + Nf ** 6 * dbd(X) * eb(X) * H(X) ** 2 * L(X) ** 3 * Qd(X)
    + -Nf ** 4 * D * ebd(X) ** 2 * H(X) ** 4 * Q(X) * Qd(X)
    + 27 * Nf ** 4 * D ** 2 * ebd(X) * H(X) ** 3 * L(X) * Q(X) * Qd(X)
    + -Nf ** 4 * B(X) * ebd(X) * H(X) ** 3 * L(X) * Q(X) * Qd(X)
    + Nf ** 4 * Bb(X) * ebd(X) * H(X) ** 3 * L(X) * Q(X) * Qd(X)
    + -Nf ** 4 * ebd(X) * G(X) * H(X) ** 3 * L(X) * Q(X) * Qd(X)
    + Nf ** 4 * ebd(X) * Gb(X) * H(X) ** 3 * L(X) * Q(X) * Qd(X)
    + -Nf ** 4 * ebd(X) * H(X) ** 4 * Hd(X) * L(X) * Q(X) * Qd(X)
    - Nf ** 3 * D ** 3 * H(X) ** 2 * L(X) ** 2 * Q(X) * Qd(X)
    + -57 * Nf ** 4 * D ** 3 * H(X) ** 2 * L(X) ** 2 * Q(X) * Qd(X)
    - (Nf ** 3 * D * B(X) * H(X) ** 2 * L(X) ** 2 * Q(X) * Qd(X)) / 2
    + -(27 * Nf ** 4 * D * B(X) * H(X) ** 2 * L(X) ** 2 * Q(X) * Qd(X)) / 2
    - (Nf ** 3 * D * Bb(X) * H(X) ** 2 * L(X) ** 2 * Q(X) * Qd(X)) / 2
    + -(19 * Nf ** 4 * D * Bb(X) * H(X) ** 2 * L(X) ** 2 * Q(X) * Qd(X)) / 2
    + 6 * Nf ** 6 * db(X) * dbd(X) * H(X) ** 2 * L(X) ** 2 * Q(X) * Qd(X)
    + -3 * Nf ** 6 * eb(X) * ebd(X) * H(X) ** 2 * L(X) ** 2 * Q(X) * Qd(X)
    - (Nf ** 3 * D * G(X) * H(X) ** 2 * L(X) ** 2 * Q(X) * Qd(X)) / 2
    + -(27 * Nf ** 4 * D * G(X) * H(X) ** 2 * L(X) ** 2 * Q(X) * Qd(X)) / 2
    - (Nf ** 3 * D * Gb(X) * H(X) ** 2 * L(X) ** 2 * Q(X) * Qd(X)) / 2
    + -(19 * Nf ** 4 * D * Gb(X) * H(X) ** 2 * L(X) ** 2 * Q(X) * Qd(X)) / 2
    - Nf ** 3 * D * H(X) ** 3 * Hd(X) * L(X) ** 2 * Q(X) * Qd(X)
    + -15 * Nf ** 4 * D * H(X) ** 3 * Hd(X) * L(X) ** 2 * Q(X) * Qd(X)
    - (3 * Nf ** 5 * D * eb(X) * H(X) * L(X) ** 3 * Q(X) * Qd(X)) / 2
    + -(25 * Nf ** 6 * D * eb(X) * H(X) * L(X) ** 3 * Q(X) * Qd(X)) / 2
    + 3 * Nf ** 6 * H(X) ** 2 * L(X) ** 3 * Ld(X) * Q(X) * Qd(X)
    + -6 * Nf ** 6 * db(X) * ebd(X) * H(X) ** 2 * L(X) * Q(X) ** 2 * Qd(X)
    - (3 * Nf ** 5 * D * db(X) * H(X) * L(X) ** 2 * Q(X) ** 2 * Qd(X)) / 2
    + -(75 * Nf ** 6 * D * db(X) * H(X) * L(X) ** 2 * Q(X) ** 2 * Qd(X)) / 2
    + (3 * Nf ** 3 * H(X) ** 2 * L(X) ** 2 * Q(X) ** 2 * Qd(X) ** 2) / 2
    + -(9 * Nf ** 6 * H(X) ** 2 * L(X) ** 2 * Q(X) ** 2 * Qd(X) ** 2) / 2
    + Nf ** 4 * D * dbd(X) * H(X) ** 4 * L(X) ** 2 * ub(X)
    + -2 * Nf ** 4 * D * ebd(X) * H(X) ** 4 * L(X) * Q(X) * ub(X)
    + 18 * Nf ** 4 * D ** 2 * H(X) ** 3 * L(X) ** 2 * Q(X) * ub(X)
    - -(Nf ** 3 * B(X) * H(X) ** 3 * L(X) ** 2 * Q(X) * ub(X)) / 2
    + (3 * Nf ** 4 * B(X) * H(X) ** 3 * L(X) ** 2 * Q(X) * ub(X)) / 2
    - -(Nf ** 3 * G(X) * H(X) ** 3 * L(X) ** 2 * Q(X) * ub(X)) / 2
    + (3 * Nf ** 4 * G(X) * H(X) ** 3 * L(X) ** 2 * Q(X) * ub(X)) / 2
    + -Nf ** 4 * H(X) ** 4 * Hd(X) * L(X) ** 2 * Q(X) * ub(X)
    - (Nf ** 5 * eb(X) * H(X) ** 2 * L(X) ** 3 * Q(X) * ub(X)) / 2
    + -(5 * Nf ** 6 * eb(X) * H(X) ** 2 * L(X) ** 3 * Q(X) * ub(X)) / 2
    - (Nf ** 5 * db(X) * H(X) ** 2 * L(X) ** 2 * Q(X) ** 2 * ub(X)) / 2
    + -(15 * Nf ** 6 * db(X) * H(X) ** 2 * L(X) ** 2 * Q(X) ** 2 * ub(X)) / 2
    - (Nf ** 3 * D ** 3 * db(X) * ebd(X) ** 2 * H(X) ** 2 * ubd(X)) / 2
    + -(21 * Nf ** 4 * D ** 3 * db(X) * ebd(X) ** 2 * H(X) ** 2 * ubd(X)) / 2
    + Nf ** 4 * D * B(X) * db(X) * ebd(X) ** 2 * H(X) ** 2 * ubd(X)
    - -(Nf ** 3 * D * Bb(X) * db(X) * ebd(X) ** 2 * H(X) ** 2 * ubd(X)) / 2
    + (3 * Nf ** 4 * D * Bb(X) * db(X) * ebd(X) ** 2 * H(X) ** 2 * ubd(X)) / 2
    + -Nf ** 4 * D * db(X) * ebd(X) ** 2 * G(X) * H(X) ** 2 * ubd(X)
    - (Nf ** 3 * D * db(X) * ebd(X) ** 2 * Gb(X) * H(X) ** 2 * ubd(X)) / 2
    + -(3 * Nf ** 4 * D * db(X) * ebd(X) ** 2 * Gb(X) * H(X) ** 2 * ubd(X)) / 2
    + Nf ** 4 * D * db(X) * ebd(X) ** 2 * H(X) ** 3 * Hd(X) * ubd(X)
    + -24 * Nf ** 4 * D ** 4 * db(X) * ebd(X) * H(X) * L(X) * ubd(X)
    + 20 * Nf ** 4 * D ** 2 * B(X) * db(X) * ebd(X) * H(X) * L(X) * ubd(X)
    + -Nf ** 4 * B(X) ** 2 * db(X) * ebd(X) * H(X) * L(X) * ubd(X)
    + 20 * Nf ** 4 * D ** 2 * Bb(X) * db(X) * ebd(X) * H(X) * L(X) * ubd(X)
    + -Nf ** 4 * B(X) * Bb(X) * db(X) * ebd(X) * H(X) * L(X) * ubd(X)
    + Nf ** 4 * Bb(X) ** 2 * db(X) * ebd(X) * H(X) * L(X) * ubd(X)
    + -12 * Nf ** 6 * D * db(X) ** 2 * dbd(X) * ebd(X) * H(X) * L(X) * ubd(X)
    - Nf ** 5 * D * db(X) * eb(X) * ebd(X) ** 2 * H(X) * L(X) * ubd(X)
    + -6 * Nf ** 6 * D * db(X) * eb(X) * ebd(X) ** 2 * H(X) * L(X) * ubd(X)
    + 20 * Nf ** 4 * D ** 2 * db(X) * ebd(X) * G(X) * H(X) * L(X) * ubd(X)
    + -2 * Nf ** 4 * B(X) * db(X) * ebd(X) * G(X) * H(X) * L(X) * ubd(X)
    + Nf ** 4 * Bb(X) * db(X) * ebd(X) * G(X) * H(X) * L(X) * ubd(X)
    + -3 * Nf ** 4 * db(X) * ebd(X) * G(X) ** 2 * H(X) * L(X) * ubd(X)
    + 20 * Nf ** 4 * D ** 2 * db(X) * ebd(X) * Gb(X) * H(X) * L(X) * ubd(X)
    + -Nf ** 4 * B(X) * db(X) * ebd(X) * Gb(X) * H(X) * L(X) * ubd(X)
    + 2 * Nf ** 4 * Bb(X) * db(X) * ebd(X) * Gb(X) * H(X) * L(X) * ubd(X)
    + -3 * Nf ** 4 * db(X) * ebd(X) * G(X) * Gb(X) * H(X) * L(X) * ubd(X)
    + 3 * Nf ** 4 * db(X) * ebd(X) * Gb(X) ** 2 * H(X) * L(X) * ubd(X)
    + -29 * Nf ** 4 * D ** 2 * db(X) * ebd(X) * H(X) ** 2 * Hd(X) * L(X) * ubd(X)
    + Nf ** 4 * B(X) * db(X) * ebd(X) * H(X) ** 2 * Hd(X) * L(X) * ubd(X)
    + -Nf ** 4 * Bb(X) * db(X) * ebd(X) * H(X) ** 2 * Hd(X) * L(X) * ubd(X)
    + Nf ** 4 * db(X) * ebd(X) * G(X) * H(X) ** 2 * Hd(X) * L(X) * ubd(X)
    + -Nf ** 4 * db(X) * ebd(X) * Gb(X) * H(X) ** 2 * Hd(X) * L(X) * ubd(X)
    + Nf ** 4 * db(X) * ebd(X) * H(X) ** 3 * Hd(X) ** 2 * L(X) * ubd(X)
    + -(Nf ** 3 * D ** 5 * db(X) * L(X) ** 2 * ubd(X)) / 2
    + (3 * Nf ** 4 * D ** 5 * db(X) * L(X) ** 2 * ubd(X)) / 2
    + -9 * Nf ** 4 * D ** 3 * B(X) * db(X) * L(X) ** 2 * ubd(X)
    + Nf ** 3 * D * B(X) ** 2 * db(X) * L(X) ** 2 * ubd(X)
    + -3 * Nf ** 4 * D * B(X) ** 2 * db(X) * L(X) ** 2 * ubd(X)
    - (Nf ** 3 * D ** 3 * Bb(X) * db(X) * L(X) ** 2 * ubd(X)) / 2
    + -(15 * Nf ** 4 * D ** 3 * Bb(X) * db(X) * L(X) ** 2 * ubd(X)) / 2
    + (Nf ** 3 * D * B(X) * Bb(X) * db(X) * L(X) ** 2 * ubd(X)) / 2
    + -(7 * Nf ** 4 * D * B(X) * Bb(X) * db(X) * L(X) ** 2 * ubd(X)) / 2
    + (Nf ** 3 * D * Bb(X) ** 2 * db(X) * L(X) ** 2 * ubd(X)) / 2
    + -(3 * Nf ** 4 * D * Bb(X) ** 2 * db(X) * L(X) ** 2 * ubd(X)) / 2
    + Nf ** 5 * D ** 2 * db(X) ** 2 * dbd(X) * L(X) ** 2 * ubd(X)
    + -14 * Nf ** 6 * D ** 2 * db(X) ** 2 * dbd(X) * L(X) ** 2 * ubd(X)
    + (Nf ** 5 * B(X) * db(X) ** 2 * dbd(X) * L(X) ** 2 * ubd(X)) / 2
    + -(3 * Nf ** 6 * B(X) * db(X) ** 2 * dbd(X) * L(X) ** 2 * ubd(X)) / 2
    + Nf ** 6 * Bb(X) * db(X) ** 2 * dbd(X) * L(X) ** 2 * ubd(X)
    + -Nf ** 5 * D ** 2 * db(X) * eb(X) * ebd(X) * L(X) ** 2 * ubd(X)
    + 14 * Nf ** 6 * D ** 2 * db(X) * eb(X) * ebd(X) * L(X) ** 2 * ubd(X)
    + -(Nf ** 5 * B(X) * db(X) * eb(X) * ebd(X) * L(X) ** 2 * ubd(X)) / 2
    + (3 * Nf ** 6 * B(X) * db(X) * eb(X) * ebd(X) * L(X) ** 2 * ubd(X)) / 2
    + -Nf ** 6 * Bb(X) * db(X) * eb(X) * ebd(X) * L(X) ** 2 * ubd(X)
    + 9 * Nf ** 4 * D ** 3 * db(X) * G(X) * L(X) ** 2 * ubd(X)
    + -Nf ** 3 * D * B(X) * db(X) * G(X) * L(X) ** 2 * ubd(X)
    + 6 * Nf ** 4 * D * B(X) * db(X) * G(X) * L(X) ** 2 * ubd(X)
    + -(Nf ** 3 * D * Bb(X) * db(X) * G(X) * L(X) ** 2 * ubd(X)) / 2
    + (7 * Nf ** 4 * D * Bb(X) * db(X) * G(X) * L(X) ** 2 * ubd(X)) / 2
    + -Nf ** 5 * db(X) ** 2 * dbd(X) * G(X) * L(X) ** 2 * ubd(X)
    + 3 * Nf ** 6 * db(X) ** 2 * dbd(X) * G(X) * L(X) ** 2 * ubd(X)
    + -(Nf ** 5 * db(X) * eb(X) * ebd(X) * G(X) * L(X) ** 2 * ubd(X)) / 2
    + (3 * Nf ** 6 * db(X) * eb(X) * ebd(X) * G(X) * L(X) ** 2 * ubd(X)) / 2
    + -2 * Nf ** 3 * D * db(X) * G(X) ** 2 * L(X) ** 2 * ubd(X)
    + 9 * Nf ** 4 * D * db(X) * G(X) ** 2 * L(X) ** 2 * ubd(X)
    - -(Nf ** 3 * D ** 3 * db(X) * Gb(X) * L(X) ** 2 * ubd(X)) / 2
    + (15 * Nf ** 4 * D ** 3 * db(X) * Gb(X) * L(X) ** 2 * ubd(X)) / 2
    + -(Nf ** 3 * D * B(X) * db(X) * Gb(X) * L(X) ** 2 * ubd(X)) / 2
    + (7 * Nf ** 4 * D * B(X) * db(X) * Gb(X) * L(X) ** 2 * ubd(X)) / 2
    + -Nf ** 3 * D * Bb(X) * db(X) * Gb(X) * L(X) ** 2 * ubd(X)
    + 4 * Nf ** 4 * D * Bb(X) * db(X) * Gb(X) * L(X) ** 2 * ubd(X)
    + -2 * Nf ** 6 * db(X) ** 2 * dbd(X) * Gb(X) * L(X) ** 2 * ubd(X)
    + Nf ** 6 * db(X) * eb(X) * ebd(X) * Gb(X) * L(X) ** 2 * ubd(X)
    + -(3 * Nf ** 3 * D * db(X) * G(X) * Gb(X) * L(X) ** 2 * ubd(X)) / 2
    + (21 * Nf ** 4 * D * db(X) * G(X) * Gb(X) * L(X) ** 2 * ubd(X)) / 2
    + -(3 * Nf ** 3 * D * db(X) * Gb(X) ** 2 * L(X) ** 2 * ubd(X)) / 2
    + (11 * Nf ** 4 * D * db(X) * Gb(X) ** 2 * L(X) ** 2 * ubd(X)) / 2
    + -45 * Nf ** 4 * D ** 3 * db(X) * H(X) * Hd(X) * L(X) ** 2 * ubd(X)
    + 10 * Nf ** 4 * D * B(X) * db(X) * H(X) * Hd(X) * L(X) ** 2 * ubd(X)
    + -7 * Nf ** 4 * D * Bb(X) * db(X) * H(X) * Hd(X) * L(X) ** 2 * ubd(X)
    + 2 * Nf ** 6 * db(X) ** 2 * dbd(X) * H(X) * Hd(X) * L(X) ** 2 * ubd(X)
    + -2 * Nf ** 6 * db(X) * eb(X) * ebd(X) * H(X) * Hd(X) * L(X) ** 2 * ubd(X)
    + 10 * Nf ** 4 * D * db(X) * G(X) * H(X) * Hd(X) * L(X) ** 2 * ubd(X)
    + -7 * Nf ** 4 * D * db(X) * Gb(X) * H(X) * Hd(X) * L(X) ** 2 * ubd(X)
    + 7 * Nf ** 4 * D * db(X) * H(X) ** 2 * Hd(X) ** 2 * L(X) ** 2 * ubd(X)
    + -5 * Nf ** 6 * D * db(X) * eb(X) * Hd(X) * L(X) ** 3 * ubd(X)
    + Nf ** 6 * db(X) * ebd(X) ** 2 * H(X) ** 2 * L(X) * Ld(X) * ubd(X)
    + -12 * Nf ** 6 * D * db(X) * ebd(X) * H(X) * L(X) ** 2 * Ld(X) * ubd(X)
    - (Nf ** 4 * D ** 2 * db(X) * L(X) ** 3 * Ld(X) * ubd(X)) / 3
    + -(28 * Nf ** 6 * D ** 2 * db(X) * L(X) ** 3 * Ld(X) * ubd(X)) / 3
    + Nf ** 6 * B(X) * db(X) * L(X) ** 3 * Ld(X) * ubd(X)
    + -(Nf ** 4 * Bb(X) * db(X) * L(X) ** 3 * Ld(X) * ubd(X)) / 3
    + (2 * Nf ** 6 * Bb(X) * db(X) * L(X) ** 3 * Ld(X) * ubd(X)) / 3
    + -Nf ** 6 * db(X) * G(X) * L(X) ** 3 * Ld(X) * ubd(X)
    + (Nf ** 4 * db(X) * Gb(X) * L(X) ** 3 * Ld(X) * ubd(X)) / 3
    + -(2 * Nf ** 6 * db(X) * Gb(X) * L(X) ** 3 * Ld(X) * ubd(X)) / 3
    + (Nf ** 4 * db(X) * H(X) * Hd(X) * L(X) ** 3 * Ld(X) * ubd(X)) / 3
    + -(5 * Nf ** 6 * db(X) * H(X) * Hd(X) * L(X) ** 3 * Ld(X) * ubd(X)) / 3
    - Nf ** 5 * D * db(X) ** 2 * ebd(X) ** 2 * H(X) * Q(X) * ubd(X)
    + -6 * Nf ** 6 * D * db(X) ** 2 * ebd(X) ** 2 * H(X) * Q(X) * ubd(X)
    + 28 * Nf ** 6 * D ** 2 * db(X) ** 2 * ebd(X) * L(X) * Q(X) * ubd(X)
    + -3 * Nf ** 6 * B(X) * db(X) ** 2 * ebd(X) * L(X) * Q(X) * ubd(X)
    + 2 * Nf ** 6 * Bb(X) * db(X) ** 2 * ebd(X) * L(X) * Q(X) * ubd(X)
    + -6 * Nf ** 6 * db(X) ** 2 * ebd(X) * G(X) * L(X) * Q(X) * ubd(X)
    + 4 * Nf ** 6 * db(X) ** 2 * ebd(X) * Gb(X) * L(X) * Q(X) * ubd(X)
    + -4 * Nf ** 6 * db(X) ** 2 * ebd(X) * H(X) * Hd(X) * L(X) * Q(X) * ubd(X)
    + 15 * Nf ** 6 * D * db(X) ** 2 * Hd(X) * L(X) ** 2 * Q(X) * ubd(X)
    - -Nf ** 3 * D ** 2 * ebd(X) ** 2 * H(X) ** 3 * Qd(X) * ubd(X)
    + 6 * Nf ** 4 * D ** 2 * ebd(X) ** 2 * H(X) ** 3 * Qd(X) * ubd(X)
    + -45 * Nf ** 4 * D ** 3 * ebd(X) * H(X) ** 2 * L(X) * Qd(X) * ubd(X)
    + 7 * Nf ** 4 * D * B(X) * ebd(X) * H(X) ** 2 * L(X) * Qd(X) * ubd(X)
    + -10 * Nf ** 4 * D * Bb(X) * ebd(X) * H(X) ** 2 * L(X) * Qd(X) * ubd(X)
    + 4 * Nf ** 6 * db(X) * dbd(X) * ebd(X) * H(X) ** 2 * L(X) * Qd(X) * ubd(X)
    + -Nf ** 6 * eb(X) * ebd(X) ** 2 * H(X) ** 2 * L(X) * Qd(X) * ubd(X)
    + 7 * Nf ** 4 * D * ebd(X) * G(X) * H(X) ** 2 * L(X) * Qd(X) * ubd(X)
    + -10 * Nf ** 4 * D * ebd(X) * Gb(X) * H(X) ** 2 * L(X) * Qd(X) * ubd(X)
    + 9 * Nf ** 4 * D * ebd(X) * H(X) ** 3 * Hd(X) * L(X) * Qd(X) * ubd(X)
    + -24 * Nf ** 4 * D ** 4 * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + 20 * Nf ** 4 * D ** 2 * B(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + -Nf ** 4 * B(X) ** 2 * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + 20 * Nf ** 4 * D ** 2 * Bb(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + -Nf ** 4 * B(X) * Bb(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + Nf ** 4 * Bb(X) ** 2 * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + -24 * Nf ** 6 * D * db(X) * dbd(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + 12 * Nf ** 6 * D * eb(X) * ebd(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + -20 * Nf ** 4 * D ** 2 * G(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + 2 * Nf ** 4 * B(X) * G(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + -Nf ** 4 * Bb(X) * G(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + 3 * Nf ** 4 * G(X) ** 2 * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + -20 * Nf ** 4 * D ** 2 * Gb(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + Nf ** 4 * B(X) * Gb(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + -2 * Nf ** 4 * Bb(X) * Gb(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + 3 * Nf ** 4 * G(X) * Gb(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + -3 * Nf ** 4 * Gb(X) ** 2 * H(X) * L(X) ** 2 * Qd(X) * ubd(X)
    - Nf ** 3 * D ** 2 * H(X) ** 2 * Hd(X) * L(X) ** 2 * Qd(X) * ubd(X)
    + -38 * Nf ** 4 * D ** 2 * H(X) ** 2 * Hd(X) * L(X) ** 2 * Qd(X) * ubd(X)
    - (Nf ** 3 * B(X) * H(X) ** 2 * Hd(X) * L(X) ** 2 * Qd(X) * ubd(X)) / 2
    + -(3 * Nf ** 4 * B(X) * H(X) ** 2 * Hd(X) * L(X) ** 2 * Qd(X) * ubd(X)) / 2
    + (Nf ** 3 * Bb(X) * H(X) ** 2 * Hd(X) * L(X) ** 2 * Qd(X) * ubd(X)) / 2
    + -(3 * Nf ** 4 * Bb(X) * H(X) ** 2 * Hd(X) * L(X) ** 2 * Qd(X) * ubd(X)) / 2
    - (Nf ** 3 * G(X) * H(X) ** 2 * Hd(X) * L(X) ** 2 * Qd(X) * ubd(X)) / 2
    + -(3 * Nf ** 4 * G(X) * H(X) ** 2 * Hd(X) * L(X) ** 2 * Qd(X) * ubd(X)) / 2
    + (Nf ** 3 * Gb(X) * H(X) ** 2 * Hd(X) * L(X) ** 2 * Qd(X) * ubd(X)) / 2
    + -(3 * Nf ** 4 * Gb(X) * H(X) ** 2 * Hd(X) * L(X) ** 2 * Qd(X) * ubd(X)) / 2
    + (Nf ** 3 * H(X) ** 3 * Hd(X) ** 2 * L(X) ** 2 * Qd(X) * ubd(X)) / 2
    + -(3 * Nf ** 4 * H(X) ** 3 * Hd(X) ** 2 * L(X) ** 2 * Qd(X) * ubd(X)) / 2
    - (Nf ** 4 * D ** 2 * eb(X) * L(X) ** 3 * Qd(X) * ubd(X)) / 3
    + -(28 * Nf ** 6 * D ** 2 * eb(X) * L(X) ** 3 * Qd(X) * ubd(X)) / 3
    + Nf ** 6 * B(X) * eb(X) * L(X) ** 3 * Qd(X) * ubd(X)
    + -(Nf ** 4 * Bb(X) * eb(X) * L(X) ** 3 * Qd(X) * ubd(X)) / 3
    + (2 * Nf ** 6 * Bb(X) * eb(X) * L(X) ** 3 * Qd(X) * ubd(X)) / 3
    + -Nf ** 6 * eb(X) * G(X) * L(X) ** 3 * Qd(X) * ubd(X)
    + (Nf ** 4 * eb(X) * Gb(X) * L(X) ** 3 * Qd(X) * ubd(X)) / 3
    + -(2 * Nf ** 6 * eb(X) * Gb(X) * L(X) ** 3 * Qd(X) * ubd(X)) / 3
    + (Nf ** 4 * eb(X) * H(X) * Hd(X) * L(X) ** 3 * Qd(X) * ubd(X)) / 3
    + -(5 * Nf ** 6 * eb(X) * H(X) * Hd(X) * L(X) ** 3 * Qd(X) * ubd(X)) / 3
    + Nf ** 5 * ebd(X) * H(X) ** 2 * L(X) ** 2 * Ld(X) * Qd(X) * ubd(X)
    + -3 * Nf ** 6 * ebd(X) * H(X) ** 2 * L(X) ** 2 * Ld(X) * Qd(X) * ubd(X)
    + Nf ** 4 * D * H(X) * L(X) ** 3 * Ld(X) * Qd(X) * ubd(X)
    - -Nf ** 5 * D * H(X) * L(X) ** 3 * Ld(X) * Qd(X) * ubd(X)
    + 10 * Nf ** 6 * D * H(X) * L(X) ** 3 * Ld(X) * Qd(X) * ubd(X)
    + -2 * Nf ** 6 * db(X) * ebd(X) ** 2 * H(X) ** 2 * Q(X) * Qd(X) * ubd(X)
    + 48 * Nf ** 6 * D * db(X) * ebd(X) * H(X) * L(X) * Q(X) * Qd(X) * ubd(X)
    + -56 * Nf ** 6 * D ** 2 * db(X) * L(X) ** 2 * Q(X) * Qd(X) * ubd(X)
    + 6 * Nf ** 6 * B(X) * db(X) * L(X) ** 2 * Q(X) * Qd(X) * ubd(X)
    + -4 * Nf ** 6 * Bb(X) * db(X) * L(X) ** 2 * Q(X) * Qd(X) * ubd(X)
    + 12 * Nf ** 6 * db(X) * G(X) * L(X) ** 2 * Q(X) * Qd(X) * ubd(X)
    + -8 * Nf ** 6 * db(X) * Gb(X) * L(X) ** 2 * Q(X) * Qd(X) * ubd(X)
    + 10 * Nf ** 6 * db(X) * H(X) * Hd(X) * L(X) ** 2 * Q(X) * Qd(X) * ubd(X)
    + -Nf ** 5 * dbd(X) * H(X) ** 2 * L(X) ** 2 * Qd(X) ** 2 * ubd(X)
    + 3 * Nf ** 6 * dbd(X) * H(X) ** 2 * L(X) ** 2 * Qd(X) ** 2 * ubd(X)
    + -6 * Nf ** 6 * ebd(X) * H(X) ** 2 * L(X) * Q(X) * Qd(X) ** 2 * ubd(X)
    - Nf ** 5 * D * H(X) * L(X) ** 2 * Q(X) * Qd(X) ** 2 * ubd(X)
    + -30 * Nf ** 6 * D * H(X) * L(X) ** 2 * Q(X) * Qd(X) ** 2 * ubd(X)
    + 9 * Nf ** 4 * D ** 2 * ebd(X) * H(X) ** 3 * L(X) * ub(X) * ubd(X)
    - -(Nf ** 3 * D ** 3 * H(X) ** 2 * L(X) ** 2 * ub(X) * ubd(X)) / 2
    + (45 * Nf ** 4 * D ** 3 * H(X) ** 2 * L(X) ** 2 * ub(X) * ubd(X)) / 2
    + -5 * Nf ** 4 * D * B(X) * H(X) ** 2 * L(X) ** 2 * ub(X) * ubd(X)
    - (Nf ** 3 * D * Bb(X) * H(X) ** 2 * L(X) ** 2 * ub(X) * ubd(X)) / 2
    + -(7 * Nf ** 4 * D * Bb(X) * H(X) ** 2 * L(X) ** 2 * ub(X) * ubd(X)) / 2
    + 2 * Nf ** 6 * db(X) * dbd(X) * H(X) ** 2 * L(X) ** 2 * ub(X) * ubd(X)
    + -Nf ** 6 * eb(X) * ebd(X) * H(X) ** 2 * L(X) ** 2 * ub(X) * ubd(X)
    + 5 * Nf ** 4 * D * G(X) * H(X) ** 2 * L(X) ** 2 * ub(X) * ubd(X)
    - -(Nf ** 3 * D * Gb(X) * H(X) ** 2 * L(X) ** 2 * ub(X) * ubd(X)) / 2
    + (7 * Nf ** 4 * D * Gb(X) * H(X) ** 2 * L(X) ** 2 * ub(X) * ubd(X)) / 2
    - -(Nf ** 3 * D * H(X) ** 3 * Hd(X) * L(X) ** 2 * ub(X) * ubd(X)) / 2
    + (9 * Nf ** 4 * D * H(X) ** 3 * Hd(X) * L(X) ** 2 * ub(X) * ubd(X)) / 2
    + -5 * Nf ** 6 * D * eb(X) * H(X) * L(X) ** 3 * ub(X) * ubd(X)
    + Nf ** 6 * H(X) ** 2 * L(X) ** 3 * Ld(X) * ub(X) * ubd(X)
    + -4 * Nf ** 6 * db(X) * ebd(X) * H(X) ** 2 * L(X) * Q(X) * ub(X) * ubd(X)
    + 30 * Nf ** 6 * D * db(X) * H(X) * L(X) ** 2 * Q(X) * ub(X) * ubd(X)
    + -6 * Nf ** 6 * H(X) ** 2 * L(X) ** 2 * Q(X) * Qd(X) * ub(X) * ubd(X)
    + 2 * Nf ** 3 * D ** 2 * db(X) ** 2 * ebd(X) ** 2 * ubd(X) ** 2
    + -(Nf ** 4 * D ** 2 * db(X) ** 2 * ebd(X) ** 2 * ubd(X) ** 2) / 2
    - (Nf ** 5 * D ** 2 * db(X) ** 2 * ebd(X) ** 2 * ubd(X) ** 2) / 2
    + -7 * Nf ** 6 * D ** 2 * db(X) ** 2 * ebd(X) ** 2 * ubd(X) ** 2
    - (Nf ** 3 * B(X) * db(X) ** 2 * ebd(X) ** 2 * ubd(X) ** 2) / 2
    + -(Nf ** 6 * B(X) * db(X) ** 2 * ebd(X) ** 2 * ubd(X) ** 2) / 2
    - (Nf ** 3 * Bb(X) * db(X) ** 2 * ebd(X) ** 2 * ubd(X) ** 2) / 4
    - -(Nf ** 4 * Bb(X) * db(X) ** 2 * ebd(X) ** 2 * ubd(X) ** 2) / 4
    - (Nf ** 5 * Bb(X) * db(X) ** 2 * ebd(X) ** 2 * ubd(X) ** 2) / 4
    + -(3 * Nf ** 6 * Bb(X) * db(X) ** 2 * ebd(X) ** 2 * ubd(X) ** 2) / 4
    + Nf ** 6 * db(X) ** 2 * ebd(X) ** 2 * G(X) * ubd(X) ** 2
    - -(Nf ** 5 * db(X) ** 2 * ebd(X) ** 2 * Gb(X) * ubd(X) ** 2) / 2
    + (3 * Nf ** 6 * db(X) ** 2 * ebd(X) ** 2 * Gb(X) * ubd(X) ** 2) / 2
    + -(Nf ** 3 * db(X) ** 2 * ebd(X) ** 2 * H(X) * Hd(X) * ubd(X) ** 2) / 2
    + (Nf ** 6 * db(X) ** 2 * ebd(X) ** 2 * H(X) * Hd(X) * ubd(X) ** 2) / 2
    + -6 * Nf ** 6 * D * db(X) ** 2 * ebd(X) * Hd(X) * L(X) * ubd(X) ** 2
    + (Nf ** 3 * db(X) ** 2 * Hd(X) ** 2 * L(X) ** 2 * ubd(X) ** 2) / 2
    + -(Nf ** 6 * db(X) ** 2 * Hd(X) ** 2 * L(X) ** 2 * ubd(X) ** 2) / 2
    - (3 * Nf ** 5 * D * db(X) * ebd(X) ** 2 * H(X) * Qd(X) * ubd(X) ** 2) / 2
    + -(15 * Nf ** 6 * D * db(X) * ebd(X) ** 2 * H(X) * Qd(X) * ubd(X) ** 2) / 2
    + 28 * Nf ** 6 * D ** 2 * db(X) * ebd(X) * L(X) * Qd(X) * ubd(X) ** 2
    + -2 * Nf ** 6 * B(X) * db(X) * ebd(X) * L(X) * Qd(X) * ubd(X) ** 2
    + 3 * Nf ** 6 * Bb(X) * db(X) * ebd(X) * L(X) * Qd(X) * ubd(X) ** 2
    + -4 * Nf ** 6 * db(X) * ebd(X) * G(X) * L(X) * Qd(X) * ubd(X) ** 2
    + 6 * Nf ** 6 * db(X) * ebd(X) * Gb(X) * L(X) * Qd(X) * ubd(X) ** 2
    + -4 * Nf ** 6 * db(X) * ebd(X) * H(X) * Hd(X) * L(X) * Qd(X) * ubd(X) ** 2
    + 12 * Nf ** 6 * D * db(X) * Hd(X) * L(X) ** 2 * Qd(X) * ubd(X) ** 2
    + -(3 * Nf ** 3 * ebd(X) ** 2 * H(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2) / 4
    + (Nf ** 4 * ebd(X) ** 2 * H(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2) / 4
    - -(Nf ** 5 * ebd(X) ** 2 * H(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2) / 4
    + (5 * Nf ** 6 * ebd(X) ** 2 * H(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2) / 4
    + -15 * Nf ** 6 * D * ebd(X) * H(X) * L(X) * Qd(X) ** 2 * ubd(X) ** 2
    + 4 * Nf ** 3 * D ** 2 * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2
    + -14 * Nf ** 6 * D ** 2 * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2
    - Nf ** 3 * B(X) * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2
    + -Nf ** 6 * B(X) * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2
    - (Nf ** 3 * Bb(X) * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2) / 2
    + -(3 * Nf ** 6 * Bb(X) * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2) / 2
    + 2 * Nf ** 6 * G(X) * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2
    + -3 * Nf ** 6 * Gb(X) * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2
    + (Nf ** 3 * H(X) * Hd(X) * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2) / 2
    + -(Nf ** 4 * H(X) * Hd(X) * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2) / 2
    + (Nf ** 5 * H(X) * Hd(X) * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2) / 2
    + -(5 * Nf ** 6 * H(X) * Hd(X) * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2) / 2
    + 12 * Nf ** 6 * D * db(X) * ebd(X) * H(X) * L(X) * ub(X) * ubd(X) ** 2
    + -Nf ** 5 * D ** 2 * db(X) * L(X) ** 2 * ub(X) * ubd(X) ** 2
    + 14 * Nf ** 6 * D ** 2 * db(X) * L(X) ** 2 * ub(X) * ubd(X) ** 2
    + -(Nf ** 5 * B(X) * db(X) * L(X) ** 2 * ub(X) * ubd(X) ** 2) / 2
    + (3 * Nf ** 6 * B(X) * db(X) * L(X) ** 2 * ub(X) * ubd(X) ** 2) / 2
    + -Nf ** 6 * Bb(X) * db(X) * L(X) ** 2 * ub(X) * ubd(X) ** 2
    + Nf ** 5 * db(X) * G(X) * L(X) ** 2 * ub(X) * ubd(X) ** 2
    + -3 * Nf ** 6 * db(X) * G(X) * L(X) ** 2 * ub(X) * ubd(X) ** 2
    + 2 * Nf ** 6 * db(X) * Gb(X) * L(X) ** 2 * ub(X) * ubd(X) ** 2
    + -2 * Nf ** 6 * db(X) * H(X) * Hd(X) * L(X) ** 2 * ub(X) * ubd(X) ** 2
    + 2 * Nf ** 6 * ebd(X) * H(X) ** 2 * L(X) * Qd(X) * ub(X) * ubd(X) ** 2
    + -12 * Nf ** 6 * D * H(X) * L(X) ** 2 * Qd(X) * ub(X) * ubd(X) ** 2
    + (Nf ** 3 * H(X) ** 2 * L(X) ** 2 * ub(X) ** 2 * ubd(X) ** 2) / 2
    + -(Nf ** 6 * H(X) ** 2 * L(X) ** 2 * ub(X) ** 2 * ubd(X) ** 2) / 2
    - (Nf * D ** 2 * ebd(X) ** 2 * H(X) ** 4 * W(X)) / 2
    + -(3 * Nf ** 2 * D ** 2 * ebd(X) ** 2 * H(X) ** 4 * W(X)) / 2
    + 18 * Nf ** 2 * D ** 3 * ebd(X) * H(X) ** 3 * L(X) * W(X)
    + -5 * Nf ** 2 * D * B(X) * ebd(X) * H(X) ** 3 * L(X) * W(X)
    + 3 * Nf ** 2 * D * Bb(X) * ebd(X) * H(X) ** 3 * L(X) * W(X)
    + -Nf ** 4 * db(X) * dbd(X) * ebd(X) * H(X) ** 3 * L(X) * W(X)
    + (Nf ** 3 * eb(X) * ebd(X) ** 2 * H(X) ** 3 * L(X) * W(X)) / 2
    + -(Nf ** 4 * eb(X) * ebd(X) ** 2 * H(X) ** 3 * L(X) * W(X)) / 2
    + 4 * Nf ** 2 * D * ebd(X) * H(X) ** 4 * Hd(X) * L(X) * W(X)
    + -(Nf * D ** 4 * H(X) ** 2 * L(X) ** 2 * W(X)) / 2
    + (37 * Nf ** 2 * D ** 4 * H(X) ** 2 * L(X) ** 2 * W(X)) / 2
    - Nf * D ** 2 * B(X) * H(X) ** 2 * L(X) ** 2 * W(X)
    + -19 * Nf ** 2 * D ** 2 * B(X) * H(X) ** 2 * L(X) ** 2 * W(X)
    + 2 * Nf ** 2 * B(X) ** 2 * H(X) ** 2 * L(X) ** 2 * W(X)
    - -(Nf * D ** 2 * Bb(X) * H(X) ** 2 * L(X) ** 2 * W(X)) / 2
    + (19 * Nf ** 2 * D ** 2 * Bb(X) * H(X) ** 2 * L(X) ** 2 * W(X)) / 2
    + -Nf ** 2 * Bb(X) ** 2 * H(X) ** 2 * L(X) ** 2 * W(X)
    - (Nf ** 3 * D * db(X) * dbd(X) * H(X) ** 2 * L(X) ** 2 * W(X)) / 2
    + -(17 * Nf ** 4 * D * db(X) * dbd(X) * H(X) ** 2 * L(X) ** 2 * W(X)) / 2
    - (Nf ** 3 * D * eb(X) * ebd(X) * H(X) ** 2 * L(X) ** 2 * W(X)) / 2
    + -(17 * Nf ** 4 * D * eb(X) * ebd(X) * H(X) ** 2 * L(X) ** 2 * W(X)) / 2
    + 2 * Nf ** 2 * G(X) ** 2 * H(X) ** 2 * L(X) ** 2 * W(X)
    + -Nf ** 2 * Gb(X) ** 2 * H(X) ** 2 * L(X) ** 2 * W(X)
    - 2 * Nf * D ** 2 * H(X) ** 3 * Hd(X) * L(X) ** 2 * W(X)
    + -24 * Nf ** 2 * D ** 2 * H(X) ** 3 * Hd(X) * L(X) ** 2 * W(X)
    + 3 * Nf ** 2 * B(X) * H(X) ** 3 * Hd(X) * L(X) ** 2 * W(X)
    - -(Nf * H(X) ** 4 * Hd(X) ** 2 * L(X) ** 2 * W(X)) / 2
    + (3 * Nf ** 2 * H(X) ** 4 * Hd(X) ** 2 * L(X) ** 2 * W(X)) / 2
    - -(3 * Nf ** 3 * D ** 2 * eb(X) * H(X) * L(X) ** 3 * W(X)) / 2
    + (33 * Nf ** 4 * D ** 2 * eb(X) * H(X) * L(X) ** 3 * W(X)) / 2
    - -Nf ** 3 * B(X) * eb(X) * H(X) * L(X) ** 3 * W(X)
    + 3 * Nf ** 4 * B(X) * eb(X) * H(X) * L(X) ** 3 * W(X)
    - Nf ** 3 * eb(X) * H(X) ** 2 * Hd(X) * L(X) ** 3 * W(X)
    + -3 * Nf ** 4 * eb(X) * H(X) ** 2 * Hd(X) * L(X) ** 3 * W(X)
    - (Nf ** 2 * eb(X) ** 2 * L(X) ** 4 * W(X)) / 8
    + (Nf ** 3 * eb(X) ** 2 * L(X) ** 4 * W(X)) / 16
    + -(Nf ** 4 * eb(X) ** 2 * L(X) ** 4 * W(X)) / 16
    - (9 * Nf ** 5 * eb(X) ** 2 * L(X) ** 4 * W(X)) / 16
    + (9 * Nf ** 6 * eb(X) ** 2 * L(X) ** 4 * W(X)) / 16
    - -(Nf ** 3 * ebd(X) * H(X) ** 3 * L(X) ** 2 * Ld(X) * W(X)) / 2
    + (3 * Nf ** 4 * ebd(X) * H(X) ** 3 * L(X) ** 2 * Ld(X) * W(X)) / 2
    - -(3 * Nf ** 3 * D * H(X) ** 2 * L(X) ** 3 * Ld(X) * W(X)) / 2
    + (17 * Nf ** 4 * D * H(X) ** 2 * L(X) ** 3 * Ld(X) * W(X)) / 2
    + -(Nf ** 3 * db(X) * ebd(X) ** 2 * H(X) ** 3 * Q(X) * W(X)) / 2
    + (Nf ** 4 * db(X) * ebd(X) ** 2 * H(X) ** 3 * Q(X) * W(X)) / 2
    + -17 * Nf ** 4 * D * db(X) * ebd(X) * H(X) ** 2 * L(X) * Q(X) * W(X)
    - (3 * Nf ** 3 * D ** 2 * db(X) * H(X) * L(X) ** 2 * Q(X) * W(X)) / 2
    + -(99 * Nf ** 4 * D ** 2 * db(X) * H(X) * L(X) ** 2 * Q(X) * W(X)) / 2
    - Nf ** 3 * B(X) * db(X) * H(X) * L(X) ** 2 * Q(X) * W(X)
    + -9 * Nf ** 4 * B(X) * db(X) * H(X) * L(X) ** 2 * Q(X) * W(X)
    - Nf ** 3 * db(X) * G(X) * H(X) * L(X) ** 2 * Q(X) * W(X)
    + -9 * Nf ** 4 * db(X) * G(X) * H(X) * L(X) ** 2 * Q(X) * W(X)
    - Nf ** 3 * db(X) * H(X) ** 2 * Hd(X) * L(X) ** 2 * Q(X) * W(X)
    + -9 * Nf ** 4 * db(X) * H(X) ** 2 * Hd(X) * L(X) ** 2 * Q(X) * W(X)
    - (3 * Nf ** 5 * db(X) * eb(X) * L(X) ** 3 * Q(X) * W(X)) / 2
    + -(9 * Nf ** 6 * db(X) * eb(X) * L(X) ** 3 * Q(X) * W(X)) / 2
    + (3 * Nf ** 3 * db(X) ** 2 * L(X) ** 2 * Q(X) ** 2 * W(X)) / 4
    + -(Nf ** 4 * db(X) ** 2 * L(X) ** 2 * Q(X) ** 2 * W(X)) / 4
    - (3 * Nf ** 5 * db(X) ** 2 * L(X) ** 2 * Q(X) ** 2 * W(X)) / 4
    + -(27 * Nf ** 6 * db(X) ** 2 * L(X) ** 2 * Q(X) ** 2 * W(X)) / 4
    - (Nf ** 3 * dbd(X) * H(X) ** 3 * L(X) ** 2 * Qd(X) * W(X)) / 2
    + -(3 * Nf ** 4 * dbd(X) * H(X) ** 3 * L(X) ** 2 * Qd(X) * W(X)) / 2
    + 3 * Nf ** 4 * ebd(X) * H(X) ** 3 * L(X) * Q(X) * Qd(X) * W(X)
    - -(3 * Nf ** 3 * D * H(X) ** 2 * L(X) ** 2 * Q(X) * Qd(X) * W(X)) / 2
    + (51 * Nf ** 4 * D * H(X) ** 2 * L(X) ** 2 * Q(X) * Qd(X) * W(X)) / 2
    - -(Nf ** 3 * H(X) ** 3 * L(X) ** 2 * Q(X) * ub(X) * W(X)) / 2
    + (9 * Nf ** 4 * H(X) ** 3 * L(X) ** 2 * Q(X) * ub(X) * W(X)) / 2
    - -(Nf ** 3 * D * db(X) * ebd(X) ** 2 * H(X) ** 2 * ubd(X) * W(X)) / 2
    + (5 * Nf ** 4 * D * db(X) * ebd(X) ** 2 * H(X) ** 2 * ubd(X) * W(X)) / 2
    + -20 * Nf ** 4 * D ** 2 * db(X) * ebd(X) * H(X) * L(X) * ubd(X) * W(X)
    + 2 * Nf ** 4 * B(X) * db(X) * ebd(X) * H(X) * L(X) * ubd(X) * W(X)
    + -Nf ** 4 * Bb(X) * db(X) * ebd(X) * H(X) * L(X) * ubd(X) * W(X)
    + 2 * Nf ** 4 * db(X) * ebd(X) * G(X) * H(X) * L(X) * ubd(X) * W(X)
    + -Nf ** 4 * db(X) * ebd(X) * Gb(X) * H(X) * L(X) * ubd(X) * W(X)
    + 2 * Nf ** 4 * db(X) * ebd(X) * H(X) ** 2 * Hd(X) * L(X) * ubd(X) * W(X)
    + -9 * Nf ** 4 * D ** 3 * db(X) * L(X) ** 2 * ubd(X) * W(X)
    - Nf ** 3 * D * B(X) * db(X) * L(X) ** 2 * ubd(X) * W(X)
    + -6 * Nf ** 4 * D * B(X) * db(X) * L(X) ** 2 * ubd(X) * W(X)
    - (Nf ** 3 * D * Bb(X) * db(X) * L(X) ** 2 * ubd(X) * W(X)) / 2
    + -(7 * Nf ** 4 * D * Bb(X) * db(X) * L(X) ** 2 * ubd(X) * W(X)) / 2
    - (Nf ** 5 * db(X) ** 2 * dbd(X) * L(X) ** 2 * ubd(X) * W(X)) / 2
    + -(3 * Nf ** 6 * db(X) ** 2 * dbd(X) * L(X) ** 2 * ubd(X) * W(X)) / 2
    - (Nf ** 5 * db(X) * eb(X) * ebd(X) * L(X) ** 2 * ubd(X) * W(X)) / 2
    + -(3 * Nf ** 6 * db(X) * eb(X) * ebd(X) * L(X) ** 2 * ubd(X) * W(X)) / 2
    - Nf ** 3 * D * db(X) * G(X) * L(X) ** 2 * ubd(X) * W(X)
    + -6 * Nf ** 4 * D * db(X) * G(X) * L(X) ** 2 * ubd(X) * W(X)
    - (Nf ** 3 * D * db(X) * Gb(X) * L(X) ** 2 * ubd(X) * W(X)) / 2
    + -(7 * Nf ** 4 * D * db(X) * Gb(X) * L(X) ** 2 * ubd(X) * W(X)) / 2
    - Nf ** 3 * D * db(X) * H(X) * Hd(X) * L(X) ** 2 * ubd(X) * W(X)
    + -15 * Nf ** 4 * D * db(X) * H(X) * Hd(X) * L(X) ** 2 * ubd(X) * W(X)
    - (Nf ** 5 * db(X) * L(X) ** 3 * Ld(X) * ubd(X) * W(X)) / 2
    + -(3 * Nf ** 6 * db(X) * L(X) ** 3 * Ld(X) * ubd(X) * W(X)) / 2
    + 3 * Nf ** 6 * db(X) ** 2 * ebd(X) * L(X) * Q(X) * ubd(X) * W(X)
    + -12 * Nf ** 4 * D * ebd(X) * H(X) ** 2 * L(X) * Qd(X) * ubd(X) * W(X)
    - Nf ** 3 * D ** 2 * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * W(X)
    + -30 * Nf ** 4 * D ** 2 * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * W(X)
    + 3 * Nf ** 4 * B(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * W(X)
    - -(Nf ** 3 * Bb(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * W(X)) / 2
    + (3 * Nf ** 4 * Bb(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * W(X)) / 2
    + -3 * Nf ** 4 * G(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * W(X)
    - (Nf ** 3 * Gb(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * W(X)) / 2
    + -(3 * Nf ** 4 * Gb(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * W(X)) / 2
    - Nf ** 3 * H(X) ** 2 * Hd(X) * L(X) ** 2 * Qd(X) * ubd(X) * W(X)
    + -3 * Nf ** 4 * H(X) ** 2 * Hd(X) * L(X) ** 2 * Qd(X) * ubd(X) * W(X)
    - (Nf ** 5 * eb(X) * L(X) ** 3 * Qd(X) * ubd(X) * W(X)) / 2
    + -(3 * Nf ** 6 * eb(X) * L(X) ** 3 * Qd(X) * ubd(X) * W(X)) / 2
    - Nf ** 5 * db(X) * L(X) ** 2 * Q(X) * Qd(X) * ubd(X) * W(X)
    + -9 * Nf ** 6 * db(X) * L(X) ** 2 * Q(X) * Qd(X) * ubd(X) * W(X)
    + Nf ** 4 * ebd(X) * H(X) ** 3 * L(X) * ub(X) * ubd(X) * W(X)
    - -(Nf ** 3 * D * H(X) ** 2 * L(X) ** 2 * ub(X) * ubd(X) * W(X)) / 2
    + (17 * Nf ** 4 * D * H(X) ** 2 * L(X) ** 2 * ub(X) * ubd(X) * W(X)) / 2
    + -2 * Nf ** 6 * db(X) * ebd(X) * L(X) * Qd(X) * ubd(X) ** 2 * W(X)
    + (Nf ** 3 * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2 * W(X)) / 2
    + -(Nf ** 4 * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2 * W(X)) / 2
    - (Nf ** 5 * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2 * W(X)) / 2
    + -(3 * Nf ** 6 * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2 * W(X)) / 2
    - (Nf ** 5 * db(X) * L(X) ** 2 * ub(X) * ubd(X) ** 2 * W(X)) / 2
    + -(3 * Nf ** 6 * db(X) * L(X) ** 2 * ub(X) * ubd(X) ** 2 * W(X)) / 2
    + (Nf * ebd(X) ** 2 * H(X) ** 4 * W(X) ** 2) / 2
    + -(Nf ** 2 * ebd(X) ** 2 * H(X) ** 4 * W(X) ** 2) / 2
    + 5 * Nf ** 2 * D * ebd(X) * H(X) ** 3 * L(X) * W(X) ** 2
    + 4 * Nf * D ** 2 * H(X) ** 2 * L(X) ** 2 * W(X) ** 2
    + -19 * Nf ** 2 * D ** 2 * H(X) ** 2 * L(X) ** 2 * W(X) ** 2
    - 2 * Nf * B(X) * H(X) ** 2 * L(X) ** 2 * W(X) ** 2
    + 4 * Nf ** 2 * B(X) * H(X) ** 2 * L(X) ** 2 * W(X) ** 2
    + -(Nf * H(X) ** 3 * Hd(X) * L(X) ** 2 * W(X) ** 2) / 2
    + (7 * Nf ** 2 * H(X) ** 3 * Hd(X) * L(X) ** 2 * W(X) ** 2) / 2
    - -Nf ** 3 * eb(X) * H(X) * L(X) ** 3 * W(X) ** 2
    + 3 * Nf ** 4 * eb(X) * H(X) * L(X) ** 3 * W(X) ** 2
    - Nf ** 3 * db(X) * H(X) * L(X) ** 2 * Q(X) * W(X) ** 2
    + -9 * Nf ** 4 * db(X) * H(X) * L(X) ** 2 * Q(X) * W(X) ** 2
    + 2 * Nf ** 4 * db(X) * ebd(X) * H(X) * L(X) * ubd(X) * W(X) ** 2
    + -Nf ** 3 * D * db(X) * L(X) ** 2 * ubd(X) * W(X) ** 2
    + 6 * Nf ** 4 * D * db(X) * L(X) ** 2 * ubd(X) * W(X) ** 2
    + -3 * Nf ** 4 * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * W(X) ** 2
    + 3 * Nf ** 2 * H(X) ** 2 * L(X) ** 2 * W(X) ** 3
    - (Nf * D ** 2 * ebd(X) ** 2 * H(X) ** 4 * Wb(X)) / 2
    + -(5 * Nf ** 2 * D ** 2 * ebd(X) ** 2 * H(X) ** 4 * Wb(X)) / 2
    + 18 * Nf ** 2 * D ** 3 * ebd(X) * H(X) ** 3 * L(X) * Wb(X)
    + -3 * Nf ** 2 * D * B(X) * ebd(X) * H(X) ** 3 * L(X) * Wb(X)
    + 5 * Nf ** 2 * D * Bb(X) * ebd(X) * H(X) ** 3 * L(X) * Wb(X)
    + -Nf ** 4 * db(X) * dbd(X) * ebd(X) * H(X) ** 3 * L(X) * Wb(X)
    - (Nf ** 3 * eb(X) * ebd(X) ** 2 * H(X) ** 3 * L(X) * Wb(X)) / 2
    + -(Nf ** 4 * eb(X) * ebd(X) ** 2 * H(X) ** 3 * L(X) * Wb(X)) / 2
    + 4 * Nf ** 2 * D * ebd(X) * H(X) ** 4 * Hd(X) * L(X) * Wb(X)
    + -(Nf * D ** 4 * H(X) ** 2 * L(X) ** 2 * Wb(X)) / 2
    + (27 * Nf ** 2 * D ** 4 * H(X) ** 2 * L(X) ** 2 * Wb(X)) / 2
    - -(Nf * D ** 2 * B(X) * H(X) ** 2 * L(X) ** 2 * Wb(X)) / 2
    + (19 * Nf ** 2 * D ** 2 * B(X) * H(X) ** 2 * L(X) ** 2 * Wb(X)) / 2
    - -(Nf * D ** 2 * Bb(X) * H(X) ** 2 * L(X) ** 2 * Wb(X)) / 2
    + (25 * Nf ** 2 * D ** 2 * Bb(X) * H(X) ** 2 * L(X) ** 2 * Wb(X)) / 2
    + -Nf ** 2 * B(X) * Bb(X) * H(X) ** 2 * L(X) ** 2 * Wb(X)
    + 6 * Nf ** 4 * D * db(X) * dbd(X) * H(X) ** 2 * L(X) ** 2 * Wb(X)
    + -6 * Nf ** 4 * D * eb(X) * ebd(X) * H(X) ** 2 * L(X) ** 2 * Wb(X)
    + Nf ** 2 * G(X) * Gb(X) * H(X) ** 2 * L(X) ** 2 * Wb(X)
    - -2 * Nf * D ** 2 * H(X) ** 3 * Hd(X) * L(X) ** 2 * Wb(X)
    + 15 * Nf ** 2 * D ** 2 * H(X) ** 3 * Hd(X) * L(X) ** 2 * Wb(X)
    + -(Nf * Bb(X) * H(X) ** 3 * Hd(X) * L(X) ** 2 * Wb(X)) / 2
    + (3 * Nf ** 2 * Bb(X) * H(X) ** 3 * Hd(X) * L(X) ** 2 * Wb(X)) / 2
    - -Nf ** 3 * D ** 2 * eb(X) * H(X) * L(X) ** 3 * Wb(X)
    + 8 * Nf ** 4 * D ** 2 * eb(X) * H(X) * L(X) ** 3 * Wb(X)
    + Nf ** 4 * Bb(X) * eb(X) * H(X) * L(X) ** 3 * Wb(X)
    + -(Nf ** 3 * ebd(X) * H(X) ** 3 * L(X) ** 2 * Ld(X) * Wb(X)) / 2
    + (3 * Nf ** 4 * ebd(X) * H(X) ** 3 * L(X) ** 2 * Ld(X) * Wb(X)) / 2
    - -Nf ** 3 * D * H(X) ** 2 * L(X) ** 3 * Ld(X) * Wb(X)
    + 6 * Nf ** 4 * D * H(X) ** 2 * L(X) ** 3 * Ld(X) * Wb(X)
    - -(Nf ** 3 * db(X) * ebd(X) ** 2 * H(X) ** 3 * Q(X) * Wb(X)) / 2
    + (Nf ** 4 * db(X) * ebd(X) ** 2 * H(X) ** 3 * Q(X) * Wb(X)) / 2
    + -12 * Nf ** 4 * D * db(X) * ebd(X) * H(X) ** 2 * L(X) * Q(X) * Wb(X)
    - Nf ** 3 * D ** 2 * db(X) * H(X) * L(X) ** 2 * Q(X) * Wb(X)
    + -24 * Nf ** 4 * D ** 2 * db(X) * H(X) * L(X) ** 2 * Q(X) * Wb(X)
    + 3 * Nf ** 4 * Bb(X) * db(X) * H(X) * L(X) ** 2 * Q(X) * Wb(X)
    + -3 * Nf ** 4 * db(X) * Gb(X) * H(X) * L(X) ** 2 * Q(X) * Wb(X)
    + (Nf ** 3 * dbd(X) * H(X) ** 3 * L(X) ** 2 * Qd(X) * Wb(X)) / 2
    + -(3 * Nf ** 4 * dbd(X) * H(X) ** 3 * L(X) ** 2 * Qd(X) * Wb(X)) / 2
    + 3 * Nf ** 4 * ebd(X) * H(X) ** 3 * L(X) * Q(X) * Qd(X) * Wb(X)
    - -Nf ** 3 * D * H(X) ** 2 * L(X) ** 2 * Q(X) * Qd(X) * Wb(X)
    + 18 * Nf ** 4 * D * H(X) ** 2 * L(X) ** 2 * Q(X) * Qd(X) * Wb(X)
    - -(Nf ** 3 * D * db(X) * ebd(X) ** 2 * H(X) ** 2 * ubd(X) * Wb(X)) / 2
    + (7 * Nf ** 4 * D * db(X) * ebd(X) ** 2 * H(X) ** 2 * ubd(X) * Wb(X)) / 2
    + -20 * Nf ** 4 * D ** 2 * db(X) * ebd(X) * H(X) * L(X) * ubd(X) * Wb(X)
    + Nf ** 4 * B(X) * db(X) * ebd(X) * H(X) * L(X) * ubd(X) * Wb(X)
    + -2 * Nf ** 4 * Bb(X) * db(X) * ebd(X) * H(X) * L(X) * ubd(X) * Wb(X)
    + Nf ** 4 * db(X) * ebd(X) * G(X) * H(X) * L(X) * ubd(X) * Wb(X)
    + -2 * Nf ** 4 * db(X) * ebd(X) * Gb(X) * H(X) * L(X) * ubd(X) * Wb(X)
    + 2 * Nf ** 4 * db(X) * ebd(X) * H(X) ** 2 * Hd(X) * L(X) * ubd(X) * Wb(X)
    + -(Nf ** 3 * D ** 3 * db(X) * L(X) ** 2 * ubd(X) * Wb(X)) / 2
    + (15 * Nf ** 4 * D ** 3 * db(X) * L(X) ** 2 * ubd(X) * Wb(X)) / 2
    - -(Nf ** 3 * D * B(X) * db(X) * L(X) ** 2 * ubd(X) * Wb(X)) / 2
    + (7 * Nf ** 4 * D * B(X) * db(X) * L(X) ** 2 * ubd(X) * Wb(X)) / 2
    - -Nf ** 3 * D * Bb(X) * db(X) * L(X) ** 2 * ubd(X) * Wb(X)
    + 4 * Nf ** 4 * D * Bb(X) * db(X) * L(X) ** 2 * ubd(X) * Wb(X)
    + -Nf ** 6 * db(X) ** 2 * dbd(X) * L(X) ** 2 * ubd(X) * Wb(X)
    + Nf ** 6 * db(X) * eb(X) * ebd(X) * L(X) ** 2 * ubd(X) * Wb(X)
    - -(Nf ** 3 * D * db(X) * G(X) * L(X) ** 2 * ubd(X) * Wb(X)) / 2
    + (7 * Nf ** 4 * D * db(X) * G(X) * L(X) ** 2 * ubd(X) * Wb(X)) / 2
    - -Nf ** 3 * D * db(X) * Gb(X) * L(X) ** 2 * ubd(X) * Wb(X)
    + 4 * Nf ** 4 * D * db(X) * Gb(X) * L(X) ** 2 * ubd(X) * Wb(X)
    - -(Nf ** 3 * D * db(X) * H(X) * Hd(X) * L(X) ** 2 * ubd(X) * Wb(X)) / 2
    + (21 * Nf ** 4 * D * db(X) * H(X) * Hd(X) * L(X) ** 2 * ubd(X) * Wb(X)) / 2
    + -Nf ** 6 * db(X) * L(X) ** 3 * Ld(X) * ubd(X) * Wb(X)
    + 2 * Nf ** 6 * db(X) ** 2 * ebd(X) * L(X) * Q(X) * ubd(X) * Wb(X)
    - -(Nf ** 3 * ebd(X) ** 2 * H(X) ** 3 * Qd(X) * ubd(X) * Wb(X)) / 2
    + (3 * Nf ** 4 * ebd(X) ** 2 * H(X) ** 3 * Qd(X) * ubd(X) * Wb(X)) / 2
    + -17 * Nf ** 4 * D * ebd(X) * H(X) ** 2 * L(X) * Qd(X) * ubd(X) * Wb(X)
    - 2 * Nf ** 3 * D ** 2 * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * Wb(X)
    + -30 * Nf ** 4 * D ** 2 * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * Wb(X)
    - (Nf ** 3 * B(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * Wb(X)) / 2
    + -(3 * Nf ** 4 * B(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * Wb(X)) / 2
    + Nf ** 3 * Bb(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * Wb(X)
    + -3 * Nf ** 4 * Bb(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * Wb(X)
    - (Nf ** 3 * G(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * Wb(X)) / 2
    + -(3 * Nf ** 4 * G(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * Wb(X)) / 2
    + Nf ** 3 * Gb(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * Wb(X)
    + -3 * Nf ** 4 * Gb(X) * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * Wb(X)
    + Nf ** 3 * H(X) ** 2 * Hd(X) * L(X) ** 2 * Qd(X) * ubd(X) * Wb(X)
    + -3 * Nf ** 4 * H(X) ** 2 * Hd(X) * L(X) ** 2 * Qd(X) * ubd(X) * Wb(X)
    + Nf ** 6 * eb(X) * L(X) ** 3 * Qd(X) * ubd(X) * Wb(X)
    + -6 * Nf ** 6 * db(X) * L(X) ** 2 * Q(X) * Qd(X) * ubd(X) * Wb(X)
    + Nf ** 4 * ebd(X) * H(X) ** 3 * L(X) * ub(X) * ubd(X) * Wb(X)
    + -6 * Nf ** 4 * D * H(X) ** 2 * L(X) ** 2 * ub(X) * ubd(X) * Wb(X)
    + 3 * Nf ** 6 * db(X) * ebd(X) * L(X) * Qd(X) * ubd(X) ** 2 * Wb(X)
    + -(Nf ** 3 * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2 * Wb(X)) / 4
    - (Nf ** 4 * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2 * Wb(X)) / 4
    + -(3 * Nf ** 5 * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2 * Wb(X)) / 4
    + (9 * Nf ** 6 * L(X) ** 2 * Qd(X) ** 2 * ubd(X) ** 2 * Wb(X)) / 4
    + -Nf ** 6 * db(X) * L(X) ** 2 * ub(X) * ubd(X) ** 2 * Wb(X)
    + 6 * Nf ** 2 * D * ebd(X) * H(X) ** 3 * L(X) * W(X) * Wb(X)
    + -Nf * D ** 2 * H(X) ** 2 * L(X) ** 2 * W(X) * Wb(X)
    + 19 * Nf ** 2 * D ** 2 * H(X) ** 2 * L(X) ** 2 * W(X) * Wb(X)
    - Nf * Bb(X) * H(X) ** 2 * L(X) ** 2 * W(X) * Wb(X)
    + -2 * Nf ** 2 * Bb(X) * H(X) ** 2 * L(X) ** 2 * W(X) * Wb(X)
    + 2 * Nf ** 4 * db(X) * ebd(X) * H(X) * L(X) * ubd(X) * W(X) * Wb(X)
    + -7 * Nf ** 4 * D * db(X) * L(X) ** 2 * ubd(X) * W(X) * Wb(X)
    - Nf ** 3 * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * W(X) * Wb(X)
    + -3 * Nf ** 4 * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * W(X) * Wb(X)
    + (Nf * ebd(X) ** 2 * H(X) ** 4 * Wb(X) ** 2) / 2
    + -(Nf ** 2 * ebd(X) ** 2 * H(X) ** 4 * Wb(X) ** 2) / 2
    + 5 * Nf ** 2 * D * ebd(X) * H(X) ** 3 * L(X) * Wb(X) ** 2
    + -(5 * Nf * D ** 2 * H(X) ** 2 * L(X) ** 2 * Wb(X) ** 2) / 2
    + (25 * Nf ** 2 * D ** 2 * H(X) ** 2 * L(X) ** 2 * Wb(X) ** 2) / 2
    - -Nf * B(X) * H(X) ** 2 * L(X) ** 2 * Wb(X) ** 2
    + Nf ** 2 * B(X) * H(X) ** 2 * L(X) ** 2 * Wb(X) ** 2
    + Nf ** 2 * Bb(X) * H(X) ** 2 * L(X) ** 2 * Wb(X) ** 2
    + -Nf * H(X) ** 3 * Hd(X) * L(X) ** 2 * Wb(X) ** 2
    + 2 * Nf ** 2 * H(X) ** 3 * Hd(X) * L(X) ** 2 * Wb(X) ** 2
    + Nf ** 4 * eb(X) * H(X) * L(X) ** 3 * Wb(X) ** 2
    + -3 * Nf ** 4 * db(X) * H(X) * L(X) ** 2 * Q(X) * Wb(X) ** 2
    + 2 * Nf ** 4 * db(X) * ebd(X) * H(X) * L(X) * ubd(X) * Wb(X) ** 2
    + -4 * Nf ** 4 * D * db(X) * L(X) ** 2 * ubd(X) * Wb(X) ** 2
    + Nf ** 3 * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * Wb(X) ** 2
    + -3 * Nf ** 4 * H(X) * L(X) ** 2 * Qd(X) * ubd(X) * Wb(X) ** 2
    - (Nf * H(X) ** 2 * L(X) ** 2 * W(X) * Wb(X) ** 2) / 2
    + -(5 * Nf ** 2 * H(X) ** 2 * L(X) ** 2 * W(X) * Wb(X) ** 2) / 2
    + (Nf * H(X) ** 2 * L(X) ** 2 * Wb(X) ** 3) / 2
    + (Nf ** 2 * H(X) ** 2 * L(X) ** 2 * Wb(X) ** 3) / 2
)
