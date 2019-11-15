#!/usr/bin/env python3

"""Babu + Leung / de Goubea + Jenkins operator field content and numerical label."""

bl_list = {
    ("L", "L", "H", "H"): 1,
    ("L", "L", "L", "eb", "H"): 2,
    ("L", "L", "Q", "db", "H"): 3,
    ("L", "L", "Q.conj", "ub.conj", "H"): 4,
    ("L", "L", "Q", "db", "H", "H", "H.conj"): 5,
    ("L", "L", "Q.conj", "ub.conj", "H", "H.conj", "H"): 6,
    ("L", "Q", "eb.conj", "Q.conj", "H", "H", "H"): 7,
    ("L", "db", "eb.conj", "ub.conj", "H"): 8,
    ("L", "L", "L", "eb", "L", "eb"): 9,
    ("L", "L", "L", "eb", "Q", "db"): 10,
    ("L", "L", "Q", "db", "Q", "db"): 11,
    ("L", "L", "Q.conj", "ub.conj", "Q.conj", "ub.conj"): 12,
    ("L", "L", "Q.conj", "ub.conj", "L", "eb"): 13,
    ("L", "L", "Q.conj", "ub.conj", "Q", "db"): 14,
    ("L", "L", "L", "db", "L.conj", "ub.conj"): 15,
    ("L", "L", "eb", "db", "eb.conj", "ub.conj"): 16,
    ("L", "L", "db", "db", "db.conj", "ub.conj"): 17,
    ("L", "L", "db", "ub", "ub.conj", "ub.conj"): 18,
    ("L", "Q", "db", "db", "eb.conj", "ub.conj"): 19,
    ("L", "db", "Q.conj", "ub.conj", "eb.conj", "ub.conj"): 20,
    ("L", "L", "L", "eb", "Q", "ub", "H", "H"): 21,
    ("L", "L", "L", "eb", "L.conj", "eb.conj", "H", "H"): 22,
    ("L", "L", "L", "eb", "Q.conj", "db.conj", "H", "H"): 23,
    ("L", "L", "Q", "db", "Q", "db", "H", "H.conj"): 24,
    ("L", "L", "Q", "db", "Q", "ub", "H", "H"): 25,
    ("L", "L", "Q", "db", "L.conj", "eb.conj", "H", "H"): 26,
    ("L", "L", "Q", "db", "Q.conj", "db.conj", "H", "H"): 27,
    ("L", "L", "Q", "db", "Q.conj", "ub.conj", "H", "H.conj"): 28,
    ("L", "L", "Q", "ub", "Q.conj", "ub.conj", "H", "H"): 29,
    ("L", "L", "L.conj", "eb.conj", "Q.conj", "ub.conj", "H", "H"): 30,
    ("L", "L", "Q.conj", "db.conj", "Q.conj", "ub.conj", "H", "H"): 31,
    ("L", "L", "Q.conj", "ub.conj", "Q.conj", "ub.conj", "H", "H.conj"): 32,
    ("eb.conj", "eb.conj", "L", "L", "eb", "eb", "H", "H"): 33,
    ("eb.conj", "eb.conj", "L", "Q", "eb", "db", "H", "H"): 34,
    ("eb.conj", "eb.conj", "L", "eb", "Q.conj", "ub.conj", "H", "H"): 35,
    ("eb.conj", "eb.conj", "Q", "db", "Q", "db", "H", "H"): 36,
    ("eb.conj", "eb.conj", "Q", "db", "Q.conj", "ub.conj", "H", "H"): 37,
    ("eb.conj", "eb.conj", "Q.conj", "ub.conj", "Q.conj", "ub.conj", "H", "H"): 38,
    ("L", "L", "L", "L", "L.conj", "L.conj", "H", "H"): 39,
    ("L", "L", "L", "Q", "L.conj", "Q.conj", "H", "H"): 40,
    ("L", "L", "L", "db", "L.conj", "db.conj", "H", "H"): 41,
    ("L", "L", "L", "ub", "L.conj", "ub.conj", "H", "H"): 42,
    ("L", "L", "L", "db", "L.conj", "ub.conj", "H", "H.conj"): 43,
    ("L", "L", "Q", "eb", "Q.conj", "eb.conj", "H", "H"): 44,
    ("L", "L", "eb", "db", "eb.conj", "db.conj", "H", "H"): 45,
    ("L", "L", "eb", "ub", "eb.conj", "ub.conj", "H", "H"): 46,
    ("L", "L", "Q", "Q", "Q.conj", "Q.conj", "H", "H"): 47,
    ("L", "L", "db", "db", "db.conj", "db.conj", "H", "H"): 48,
    ("L", "L", "db", "ub", "db.conj", "ub.conj", "H", "H"): 49,
    ("L", "L", "db", "db", "db.conj", "ub.conj", "H", "H.conj"): 50,
    ("L", "L", "ub", "ub", "ub.conj", "ub.conj", "H", "H"): 51,
    ("L", "L", "db", "ub", "ub.conj", "ub.conj", "H", "H.conj"): 52,
    ("L", "L", "db", "db", "ub.conj", "ub.conj", "H.conj", "H.conj"): 53,
    ("L", "Q", "Q", "db", "Q.conj", "eb.conj", "H", "H"): 54,
    ("L", "Q", "Q.conj", "Q.conj", "eb.conj", "ub.conj", "H", "H"): 55,
    ("L", "Q", "db", "db", "eb.conj", "db.conj", "H", "H"): 56,
    ("L", "db", "Q.conj", "ub.conj", "eb.conj", "db.conj", "H", "H"): 57,
    ("L", "ub", "Q.conj", "ub.conj", "eb.conj", "ub.conj", "H", "H"): 58,
    ("L", "Q", "db", "db", "eb.conj", "ub.conj", "H", "H.conj"): 59,
    ("L", "db", "Q.conj", "ub.conj", "eb.conj", "ub.conj", "H", "H.conj"): 60,
    ("L", "L", "H", "H", "L", "eb", "H.conj"): 61,
    ("L", "L", "L", "eb", "H", "L", "eb", "H.conj"): 62,
    ("L", "L", "Q", "db", "H", "L", "eb", "H.conj"): 63,
    ("L", "L", "Q.conj", "ub.conj", "H", "L", "eb", "H.conj"): 64,
    ("L", "eb.conj", "ub.conj", "db", "H", "L", "eb", "H.conj"): 65,
    # ("L", "L", "H", "H", "Q", "db", "H.conj"): 66,                    # like 5
    # ("L", "L", "L", "eb", "H", "Q", "db", "H.conj"): 67,              # like 63
    # ("L", "L", "Q", "db", "H", "Q", "db", "H.conj"): 68,              # like 24
    # ("L", "L", "Q.conj", "ub.conj", "H", "Q", "db", "H.conj"): 69,    # like 28
    # ("L", "eb.conj", "ub.conj", "db", "H", "Q", "db", "H.conj"): 70,  # like 59
    ("L", "L", "H", "H", "Q", "ub", "H"): 71,
    # ("L", "L", "L", "eb", "H", "Q", "ub", "H"): 72,                   # like 21
    # ("L", "L", "Q", "db", "H", "Q", "ub", "H"): 73,                   # like 25
    # ("L", "L", "Q.conj", "ub.conj", "H", "Q", "ub", "H"): 74,         # like 29
    ("L", "eb.conj", "ub.conj", "db", "H", "Q", "ub", "H"): 75,
    ("ub.conj", "ub.conj", "db", "db", "eb.conj", "eb.conj"): 76,
}

BL_LIST = {}
for k, v in bl_list.items():
    BL_LIST[tuple(sorted(k))] = v
