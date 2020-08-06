#!/usr/bin/env python

import string
from sympy.core.numbers import Zero

# Be careful, some of these (e.g. \x and \r) cause various issues
TEX_GREEK_LOWERCASE = [
    r"\alpha",
    r"\beta",
    r"\gamma",
    r"\delta",
    r"\epsilon",
    r"\zeta",
    r"\eta",
    r"\theta",
    r"\iota",
    r"\kappa",
    r"\lambda",
    r"\mu",
    r"\nu",
    r"\xi",
    r"\pi",
    r"\rho",
    r"\sigma",
    r"\tau",
    r"\upsilon",
    r"\phi",
    r"\chi",
    r"\psi",
    r"\omega",
]

DOTTED_TEX_GREEK_LOWERCASE = [rf"\dot{{{a}}}" for a in TEX_GREEK_LOWERCASE]


def is_deriv_in(coll: list):
    """Check to see if there are any derivatives acting on fields in `coll`."""
    for f in coll:
        if str(f)[0] == "D":
            return True

    return False


def is_fieldstrength_in(coll: list):
    """Check to see if there are any field strengths in `coll`."""
    first = lambda f: str(f).split("D")[-1][0]  # first non deriv
    for f in coll:
        if first(f) == "B" or first(f) == "G" or first(f) == "W":
            return True

    return False


def is_linear_in_deriv(coll: list):
    """Check to see if `coll` is linear in D."""

    if is_fieldstrength_in(coll):
        return False

    fields_string = " ".join(str(f) for f in coll)
    return fields_string.count("D") == 1


def strip_parens(s: str):
    """If a string is surrounded in parens, remove them."""
    if s[0] == "(" and s[-1] == ")":
        return s[1:-1]

    return s


def repr_tree(expr, string="", spaces=2):
    """Returns a string representation of a nested tree structure of tuples.

    """
    space = spaces * " "

    if not isinstance(expr, tuple):
        return expr.__repr__()

    node, left, right = expr

    string += node.__repr__() + "\n"
    string += space + repr_tree(right, spaces=spaces + 2) + "\n"
    string += space + repr_tree(left, spaces=spaces + 2)

    return string


def to_tex(label):
    char = str(label)[0]

    if char in string.ascii_letters:
        return label

    greek = "αβγδεζηθικλμνξπρστυφχψω"
    d = dict(zip(greek, TEX_GREEK_LOWERCASE))
    return d[char]


def safe_nocoeff(expr):
    if isinstance(expr, Zero) or isinstance(expr, int):
        return 0
    return expr.nocoeff
