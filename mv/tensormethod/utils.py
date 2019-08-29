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
