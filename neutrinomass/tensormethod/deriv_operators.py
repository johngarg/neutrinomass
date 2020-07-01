#!/usr/bin/env python3

import sympy

from typing import List, Dict
from itertools import groupby
from string import ascii_lowercase
from sympy import Function
from sympy import symbols
from sympy.abc import X

from neutrinomass.tensormethod.core import Field, Operator
from neutrinomass.tensormethod.contract import invariants
from neutrinomass.tensormethod.parse_hs import parse


D = symbols("D")

H, Hd, L, Ld, Q, Qd = symbols("H Hd L Ld Q Qd", cls=Function)
eb, ebd, ub, ubd, db, dbd = symbols("eb ebd ub ubd db dbd", cls=Function)
G, Gb, W, Wb, B, Bb = symbols("G Gb W Wb B Bb", cls=Function)

# From
# >>> npoint_fieldstrings(
#         n,
#         func=lambda c: sum(f.charges["l"] for f in c) == 2,
#         derivs=True
#     )
deriv_operator_names = {
    # LLdbDub†(00000)(0)
    "1": D * L(X) ** 2 * db(X) * ubd(X),
    # LLHHDD(00000)(0)
    "2": D ** 2 * L(X) ** 2 * H(X) ** 2,
    # LHHHDeb†(00000)(0)
    "3": D * L(X) * ebd(X) * H(X) ** 3,
    # LLLebHDD(00000)(0)
    "4": D ** 2 * L(X) ** 3 * eb(X) * H(X),
    # LLLHHDL†(00000)(0)
    "5": D * L(X) ** 3 * H(X) ** 2 * Ld(X),
    # LLebHHDeb†(00000)(0)
    "6": D * L(X) ** 2 * H(X) ** 2 * eb(X) * ebd(X),
    # LLQdbHDD(00000)(0)
    "7": D ** 2 * L(X) ** 2 * Q(X) * db(X) * H(X),
    # LLQHHDQ†(00000)(0)
    "8": D * L(X) ** 2 * Q(X) * H(X) ** 2 * Qd(X),
    # LLdbHHDdb†(00000)(0)
    "9": D * L(X) ** 2 * db(X) * H(X) ** 2 * dbd(X),
    # LLdbHDub†H†(00000)(0)
    "10": D * L(X) ** 2 * db(X) * H(X) * ubd(X) * Hd(X),
    # LLdbDDDub†(00000)(0)
    "11": D ** 3 * L(X) ** 2 * db(X) * ubd(X),
    # LLubHHDub†(00000)(0)
    "12": D * L(X) ** 2 * ub(X) * H(X) ** 2 * ubd(X),
    # LLHDDQ†ub†(00000)(0)
    "13": D ** 2 * L(X) ** 2 * H(X) * Qd(X) * ubd(X),
    # LQdbHHDeb†(00000)(0)
    "14": D * L(X) * Q(X) * db(X) * H(X) ** 2 * ebd(X),
    # LdbHDDeb†ub†(00000)(0)
    "15": D ** 2 * L(X) * db(X) * H(X) * ebd(X) * ubd(X),
    # LHHDeb†Q†ub†(00000)(0)
    "16": D * L(X) * H(X) ** 2 * ebd(X) * Qd(X) * ubd(X),
    # dbHHDeb†eb†ub†(00000)(0)
    "17": D * db(X) * H(X) ** 2 * ebd(X) ** 2 * ubd(X),
    # LLHHHDDH†(00000)(0),
    "18": D ** 2 * L(X) ** 2 * H(X) ** 3 * Hd(X),
    # LLHHDDDD(00000)(0),
    "19": D ** 4 * L(X) ** 2 * H(X) ** 2,
    # LHHHHDeb†H†(00000)(0),
    "20": D * L(X) * ebd(X) * Hd(X) * H(X) ** 4,
    # LHHHDDDeb†(00000)(0),
    "21": D ** 3 * L(X) * H(X) ** 3 * ebd(X),
    # HHHHDDeb†eb†(00000)(0)
    "22": D ** 2 * ebd(X) ** 2 * H(X) ** 4,
}

# need to remove operators with unwanted derivative structures
def remove_unwanted_deriv_structures(fieldstrings: List[Field]) -> List[Field]:
    out = []
    for fieldstring in fieldstrings:
        for field in fieldstring:
            if field.derivs and field.lorentz_irrep not in ((1, 1), (0, 1), (1, 0)):
                break

            if field.derivs == 2:
                break

        else:
            out.append(fieldstring)
            continue

    return out


def keep_longest_invariants(fieldstrings: List[Field]) -> List[Operator]:
    longest_len = 0
    longest = None
    for fieldstring in fieldstrings:
        invs = invariants(*fieldstring)
        if len(invs) > longest_len:
            longest = invs

    if longest is None:
        return []

    return longest


# returns the invariant that we will put in the paper's table of derivative
# operators
def longest_deriv_invariants(
    operator_dict: Dict[str, sympy.mul.Mul]
) -> Dict[str, List[Operator]]:
    """Take pattern of derivatives acting on fields that maximizes non-vanishing
    structures so that they can all be names.

    """
    out = {}
    for k, v in operator_dict.items():
        fieldstrings = remove_unwanted_deriv_structures(parse([v]))
        out[k] = keep_longest_invariants(fieldstrings)

    return out


def label_operators(
    operator_dict: Dict[str, List[Operator]]
) -> Dict[str, List[Operator]]:
    """Returns a flattened dictionary mapping indexed operator labels to operators.

    Example:
       >>> label_operators(longest_deriv_invariants(deriv_operator_names))

    """
    out = {}
    for k, v in operator_dict.items():
        for i, op in enumerate(v):
            label = k + ascii_lowercase[i] if len(v) > 1 else k
            out[label] = op
    return out
