#!/usr/bin/env python

"""Functions to multiply tensors by epsilons and deltas."""

from typing import List
from core import Index, eps, decompose_product
from sm import Q, H, L, eb, db, ub


def contract(left: "Field", right: "Field", out: str):
    """
    Example:
        >>> contract(left=Field("L", "10001"), right=Field("L", "10001"), out="00000")
        L(u0, i0) * L(u1, i1) * eps(-i0, -i1) * eps(-u0, -u1)

    """
    # Create an index dict mapping index type to position in dynkin str
    index_dict = {}
    for pos, (_, label) in enumerate(Index.get_dynkin_labels()):
        index_dict[label] = pos

    left_indexed, right_indexed = left.fresh_indices(), right.fresh_indices()

    prod = 1
    for index_type in ["u", "d", "i"]:
        n_target_indices = int(out[index_dict[index_type]])
        n_input_indices = (
            left.dynkin_ints[index_dict[index_type]]
            + right.dynkin_ints[index_dict[index_type]]
        )

        if n_target_indices == n_input_indices:
            continue
        elif (
            n_target_indices < n_input_indices
            and (n_target_indices - n_input_indices) % 2 == 0
        ):
            # get indices of index_type from tensors
            left_indices = left_indexed.indices_by_type[
                Index.get_index_types()[index_type]
            ]
            right_indices = right_indexed.indices_by_type[
                Index.get_index_types()[index_type]
            ]
            eps_indices = list(zip(left_indices, right_indices))

            counter = n_input_indices
            for i, j in eps_indices:
                prod *= eps(" ".join([(-i).label, (-j).label]))
                counter -= 2
                if counter == n_target_indices:
                    break
        else:
            raise Exception("Can't perform the contraction with epsilons.")

    return left_indexed * right_indexed * prod
