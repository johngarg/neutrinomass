#!/usr/bin/env python

"""Functions to multiply tensors by epsilons and deltas."""

from typing import List, Union
from core import (
    Field,
    IndexedField,
    Index,
    eps,
    delta,
    Prod,
    Operator,
    decompose_product,
)
from sm import Q, H, L, eb, db, ub
import itertools
from itertools import permutations
from kanren import run, eq, membero, var, conde

Contractable = Union[Field, IndexedField, Operator]
Indexed = Union[IndexedField, Operator]


def colour_singlets(operators: List[Operator]):
    """Contracts colour indices into singlets."""
    result = []
    for op in operators:

        # collect colour indices
        colour_indices = []
        for i in op.free_indices:
            if i.index_type == Index.get_index_types()["c"]:
                colour_indices.append(i)

        # separate up and down indices
        ups, downs = [], []
        for i in colour_indices:
            if i.is_up:
                ups.append(i)
            else:
                downs.append(i)

        # can only contract into singlet if equal number of raised and lowered
        # colour indices
        if len(ups) != len(downs):
            raise ValueError("Cannot contract colour indices into a singlet.")

        delta_index_combos = [tuple(zip(perm, downs)) for perm in permutations(ups)]
        for combo in delta_index_combos:
            deltas = [delta(" ".join([(-j).label, (-i).label])) for i, j in combo]

            # multiply deltas into operator
            prod = op
            for d in deltas:
                prod *= d
            result.append(prod)

    return result


def contract_su2_helper(
    left: Indexed, right: Indexed, out: Union[int, str], index_type: str
):
    """Returns epsilon structures so that `left * right * epsilon` for each
    structure transforms as `out`.

    Example:
        >>> contract(left=L("u0 i0"), right=L("u1 i1"), out=0, index_type="i")
        [eps(-i0, -i1)]

    """
    # create an index dict mapping index type to position in dynkin str
    index_dict = {}
    for pos, (_, label) in enumerate(Index.get_dynkin_labels()):
        index_dict[label] = pos

    # if out is a string, take to be dynkin string and extract appropriate
    # dynkin label for index type
    if isinstance(out, int):
        n_target_indices = out
    else:
        n_target_indices = int(out[index_dict[index_type]])

    n_input_indices = (
        left.dynkin_ints[index_dict[index_type]]
        + right.dynkin_ints[index_dict[index_type]]
    )

    assert n_target_indices <= n_input_indices
    assert (n_target_indices - n_input_indices) % 2 == 0

    # if target indices == input indices, no need for epsilons
    if n_target_indices == n_input_indices:
        return [1]

    # get indices of index_type from tensors
    left_indices = left.indices_by_type[Index.get_index_types()[index_type]]
    right_indices = right.indices_by_type[Index.get_index_types()[index_type]]

    # There is in general more than one way to contract the indices.
    # Construct every possible way and return a list of results.

    # get different epsilon combinations
    shortest, longest = sorted([left_indices, right_indices], key=len)
    eps_index_combos = [tuple(zip(perm, shortest)) for perm in permutations(longest)]

    results = []
    for combo in eps_index_combos:
        prod = 1
        counter = n_input_indices
        for i, j in combo:
            prod *= eps(" ".join([(-i).label, (-j).label]))
            counter -= 2
            if counter == n_target_indices:
                results.append(prod)

    return results


def contract_su2(left: Contractable, right: Contractable, out: str, ignore=[]):
    """Contracts the SU2 structure into an operator transforming like `out`.

    """
    # make indexed objects
    if not isinstance(left, (IndexedField, Operator)):
        left = left.fresh_indices()
    if not isinstance(right, (IndexedField, Operator)):
        right = right.fresh_indices()

    epsilons = []
    for i in ["u", "d", "i"]:

        # optional out to ignore contracting index structure
        if i in ignore:
            epsilons.append([1])
            continue

        junk = contract_su2_helper(left, right, out, i)
        epsilons.append(junk)

    return [left * right * u * d * i for u, d, i in itertools.product(*epsilons)]


# def contract(left: Contractable, right: Contractable, out: str):
#     """
#     Example:
#         >>> contract(left=Field("L", "10001"), right=Field("L", "10001"), out="00000")
#         L(u0, i0) * L(u1, i1) * eps(-i0, -i1) * eps(-u0, -u1)

#     For now, contracts SU(2) structures only.

#     """
#     # Create an index dict mapping index type to position in dynkin str
#     index_dict = {}
#     for pos, (_, label) in enumerate(Index.get_dynkin_labels()):
#         index_dict[label] = pos

#     if not isinstance(left, (IndexedField, Operator)):
#         left_indexed = left.fresh_indices()
#     else:
#         left_indexed = left

#     if not isinstance(right, (IndexedField, Operator)):
#         right_indexed = right.fresh_indices()
#     else:
#         right_indexed = right

#     results = [1]
#     for index_type in ["u", "d", "i"]:
#         n_target_indices = int(out[index_dict[index_type]])
#         n_input_indices = (
#             left.dynkin_ints[index_dict[index_type]]
#             + right.dynkin_ints[index_dict[index_type]]
#         )

#         if n_target_indices == n_input_indices:
#             index_results = [1]
#         elif (
#             n_target_indices < n_input_indices
#             and (n_target_indices - n_input_indices) % 2 == 0
#         ):
#             # get indices of index_type from tensors
#             left_indices = left_indexed.indices_by_type[
#                 Index.get_index_types()[index_type]
#             ]
#             right_indices = right_indexed.indices_by_type[
#                 Index.get_index_types()[index_type]
#             ]

#             # There is in general more than one way to contract the indices.
#             # Construct every possible way and return a list of results.

#             # get different epsilon combinations
#             shortest, longest = sorted([left_indices, right_indices], key=len)
#             eps_index_combos = [
#                 tuple(zip(perm, shortest)) for perm in permutations(longest)
#             ]

#             index_results = []
#             for combo in eps_index_combos:
#                 prod = 1
#                 counter = n_input_indices
#                 for i, j in combo:
#                     prod *= eps(" ".join([(-i).label, (-j).label]))
#                     counter -= 2
#                     if counter == n_target_indices:
#                         index_results.append(prod)

#         else:
#             raise Exception("Can't perform the contraction with epsilons.")

#         results = product(index_results, results)

#     breakpoint()
#     return [left_indexed * right_indexed * res for res in results]


def construct_operators(expr_tree: Prod, _stack=[], ignore=[]):
    """Constructs Operator objects from product tree. There will be more than one
    since there will be many ways to contract the indices at each step.

    """
    # Algorithm logic: unload tree into a stack and read off elements, calling
    # contract on them and building up the operators.

    while isinstance(expr_tree, Prod):
        irrep, left, right = expr_tree
        _stack.append(irrep)
        _stack.append(right)
        expr_tree = left

    results = [expr_tree]  # starts as first field
    while _stack:
        field = _stack.pop()
        irrep = _stack.pop()

        updated_results = []
        for result in results:
            updated_results += contract_su2(
                left=result, right=field, out=irrep.dynkin, ignore=ignore
            )

        results = updated_results

    return results


def invariants(*fields, ignore=[]):
    prod = decompose_product(*fields)
    singlets = [i for i in prod if i.is_singlet]

    if not singlets:
        return []

    result = []
    for s in singlets:
        operators = construct_operators(s.walked(), ignore=ignore)

        if "c" not in ignore:
            operators = colour_singlets(operators)

        result += operators

    return result


# def old_delete_duplicate_operators(operators):
#     to_remove = set()
#     counter = 0
#     for i in range(len(operators)):
#         for j in range(i + 1, len(operators)):
#             # TODO Shouldn't need to interact with sympy objects, fix the below
#             if isinstance(operators[i], Operator):
#                 a = operators[i].simplify()
#                 operators[i] = a
#             else:
#                 a = operators[i]

#             if isinstance(operators[j], Operator):
#                 b = operators[j].simplify()
#                 operators[j] = b
#             else:
#                 b = operators[j]

#             if (a == 0 and b != 0) or (b == 0 and a != 0):
#                 continue
#             elif a == 0 and b == 0:
#                 to_remove.add(j)
#                 continue
#             elif a.nocoeff == b.nocoeff:
#                 to_remove.add(j)

#     return [op for i, op in enumerate(operators) if i not in to_remove]


def compare_operators(left, right):
    if left.nocoeff == right.nocoeff:
        return True, []

    mapping = []

    # epsilons are sufficient to look at
    left_eps = [t for t in left.args if str(t).startswith("Eps")]
    right_eps = [t for t in right.args if str(t).startswith("Eps")]

    for l, r in zip(left_eps, right_eps):
        left_label, left_indices = l.args
        right_label, right_indices = r.args

        if left_label != right_label:
            return False, []

        # remove minus sign differences and other potential ambiguities by
        # sorting epsilon indices
        left_indices, right_indices = sorted(left_indices), sorted(right_indices)

        # walk down indices saving mappings that make them equal
        if left_indices != right_indices:
            for i, j in zip(left_indices, right_indices):
                if i != j:
                    mapping.append((i, j))

    return True, mapping


def is_identity_mapping(mapping, op):
    for i, j in mapping:

        # disregard free indices
        if isinstance(i, Index) and isinstance(j, Index):
            continue

        for k, v in op.duplicate_fields.items():
            # TODO May lead to problems for colour duplicates because of raising
            # and lowering
            if (-i in v) and (-j in v):
                break
        else:
            return False
    return True


def clean_operators(operators):
    """Delete duplicates (ignoring relabellings) and remove vanishing operators."""
    seen, out = set(), []
    for op in operators:
        # TODO Rewrite simplify so that it still returns a regular operator
        op_tens_mul = op.simplify()
        if not op_tens_mul:
            continue
        if not op_tens_mul in seen:
            out.append(op)
            seen.add(op_tens_mul)

    return out


def remove_relabellings(operators):
    to_remove = set()
    counter = 0
    for i in range(len(operators)):
        for j in range(i + 1, len(operators)):
            possible_match, mapping = compare_operators(operators[i], operators[j])
            if not possible_match:
                continue

            if is_identity_mapping(mapping, operators[i]):
                to_remove.add(j)

    return [op for i, op in enumerate(operators) if i not in to_remove]
