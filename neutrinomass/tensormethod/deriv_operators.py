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

# Dimension 7
# [L, eb.conj, H, H, D(H, "11")]
od1 = D * L(X) * ebd(X) * H(X) ** 3

# [L, D(L, "01"), H, D(H, "11")]
od2 = D ** 2 * L(X) ** 2 * H(X) ** 2

# Dimension 9
# [L, L, H, H, H, H.conj, "D", "D"]
od3 = D ** 2 * L(X) * L(X) * H(X) ** 3 * Hd(X)

# [L, L, L, L.conj, H, H, "D"]
od4 = D * L(X) ** 3 * Ld(X) * H(X) ** 2

# [L, eb.conj, D(H, "11"), H, H, H, H.conj]
od5 = D * L(X) * ebd(X) * H(X) ** 4 * Hd(X)

# [L, L, eb.conj, eb, H, H, "D"]
od6 = D * L(X) ** 2 * ebd(X) * eb(X) * H(X) ** 2

# [eb.conj, eb.conj, D(H, "11"), D(H, "11"), H, H]
od7 = D ** 2 * ebd(X) ** 2 * H(X) ** 4

# [L, L, H, H, Q, Q.conj, "D"]
od8 = D * Q(X) * Qd(X) * L(X) ** 2 * H(X) ** 2

# [L, L, H, H, ub, ub.conj, "D"]
od9 = D * L(X) ** 2 * H(X) ** 2 * ub(X) * ubd(X)

# [L, L, H, H, db, db.conj, "D"]
od10 = D * L(X) ** 2 * H(X) ** 2 * db(X) * dbd(X)

# [L, eb.conj, Q.conj, ub.conj, H, H, "D"]
od11 = D * L(X) * ebd(X) * Qd(X) * ubd(X) * H(X) ** 2

# [L, Q, eb.conj, db, H, H, "D"]
od12 = D * L(X) * Q(X) * ebd(X) * db(X) * H(X) ** 2

# [eb.conj, eb.conj, db, ub.conj, H, D(H, "11")]
od13 = D * ebd(X) * ebd(X) * db(X) * ubd(X) * H(X) ** 2

# [L, L, db, ub.conj, H, H.conj, "D"]
od14 = D * L(X) ** 2 * db(X) * ubd(X) * H(X) * Hd(X)

deriv_operator_names = {
    "1": od1,
    "2": od2,
    "3": od3,
    "4": od4,
    "5": od5,
    "6": od6,
    "7": od7,
    "8": od8,
    "9": od9,
    "10": od10,
    "11": od11,
    "12": od12,
    "13": od13,
    "14": od14,
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
