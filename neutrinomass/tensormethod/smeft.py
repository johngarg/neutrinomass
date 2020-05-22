#!/usr/bin/env python3

"""Functions to generate ΔL = 2 SMEFT operators"""

from typing import List
from itertools import combinations_with_replacement

from neutrinomass.tensormethod.sm import L, Q, H, eb, ub, db
from neutrinomass.tensormethod.core import Field, Operator, Prod, decompose_product

# from neutrinomass.tensormethod.contract import construct_operators


# add lepton number
L.charges["L"] = 1
Q.charges["L"] = 0
H.charges["L"] = 0
ub.charges["L"] = 0
db.charges["L"] = 0
eb.charges["L"] = -1

SM = [L, Q, H, eb, ub, db]
FIELDS = SM + [f.conj for f in SM]


def prod_mass_dim(prod: Prod) -> int:
    """Returns the mass dimension of a Prod object"""
    out = []
    while True:
        irrep, left, right = prod
        out.append(right)

        if not isinstance(left, Prod):
            out.append(left)
            return sum(map(lambda x: x.mass_dim, out))

        prod = left


def generate_operators(
    max_dim: int, fields: List[Field], verbose: bool = False
) -> List[Operator]:
    """Returns a list of operators of maximum dimension `max_dim` built out of
    `fields`.

    """
    # max number of allowed fields = max dim since lowest dim fields are scalars
    out = []
    for n_fields in range(2, max_dim + 1):
        combos = combinations_with_replacement(fields, n_fields)
        for combo in combos:
            products = decompose_product(*combo)
            for prod in products:
                if prod.is_singlet and prod_mass_dim(prod.walked()) <= max_dim:
                    if verbose:
                        print(prod)
                    out.append(prod)

    return out
