#!/usr/bin/env python


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


TEX_GREEK_LOWERCASE = [
    r"\alpha",
    r"\kappa",
    r"\upsilon",
    r"\beta",
    r"\zeta",
    r"\lambda",
    r"\pi",
    r"\phi",
    r"\gamma",
    r"\eta",
    r"\mu",
    r"\rho",
    r"\chi",
    r"\delta",
    r"\theta",
    r"\nu",
    r"\sigma",
    r"\psi",
    r"\epsilon",
    r"\iota",
    r"\xi",
    r"\tau",
    r"\omega",
]

DOTTED_TEX_GREEK_LOWERCASE = [rf"\dotted{{{a}}}" for a in TEX_GREEK_LOWERCASE]
