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
from sm import *
import itertools
from itertools import combinations, permutations, product
from functools import reduce

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


def unsimplified_invariants(*fields, ignore=[]):
    """Generates invariant operators containing ``fields`` and ignoring the
    structures in ``ignore``. The structures are not filtered, simplified or
    processed at all.

    """
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


def extract_relabellings(x, y):
    """Takes two structures of the same index type and returns a list of lists of
    tuples representing index substitutions on x that would make it equal to y.

    Example:
        >>> extract_relabellings([Eps(-I_0, -I_2), Eps(-I_3, -I_1)], [Eps(-I_3, -I_2), Eps(-I_2, -I_1)])
        [[(I_0, I_3), (I_3, I_2)],
         [(I_1, I_2), (I_0, I_1)],
         [(I_2, I_3), (I_0, I_2), (I_3, I_1), (I_1, I_2)],
         ...
        ]

    """
    perms = permutations(sorted(pair.indices) for pair in x)
    y_indices = [sorted(pair.indices) for pair in y]
    mappings = []
    for perm in perms:
        mapping = []
        for (i, j), (k, l) in zip(perm, y_indices):
            if i != k:
                mapping.append((i, k))

            if j != l:
                mapping.append((j, l))

        mappings.append(mapping)

    return mappings


def compare_operators(left, right):
    if left.nocoeff == right.nocoeff:
        return True, []

    # Epsilons should be sufficient to look at.
    #
    # Separate into lists of lists of invariant symbols:
    # e.g.
    # left_structures = [[Eps(-I_0, -I_2), Eps(-I_3, -I_1)],
    #                    [Eps(-U_0, -U_3), Eps(-U_1, -U_4), Eps(-U_5, -U_2)],
    #                    ...]
    #
    # Walk through these for left and right, but take a cartesian product when
    # comparing
    left_structures = left.by_structures()
    right_structures = right.by_structures()
    assert len(left_structures) == len(right_structures)

    relabellings_by_structures = [None] * len(left_structures)

    for i, (ls, rs) in enumerate(zip(left_structures, right_structures)):
        relabellings = extract_relabellings(ls, rs)
        relabellings_by_structures[i] = relabellings

    # Is the right thing to do to return the structures separately but only say
    # there is an identity mapping if both substitutions simultaneously give you
    # the same operator?
    return [
        reduce(lambda x, y: x + y, structures)
        for structures in product(*relabellings_by_structures)
    ]


def is_identity_mapping(mapping, op):
    """This and function below are currently broken for more than one structure. See
    comments on invariants function.

    """
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


def exists_identity_mapping(mappings, op):
    for mapping in mappings:
        if is_identity_mapping(mapping, op):
            return True
    return False


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
            mappings = compare_operators(operators[i], operators[j])

            if exists_identity_mapping(mappings, operators[i]):
                to_remove.add(j)

    return [op for i, op in enumerate(operators) if i not in to_remove]


def invariants(*fields, ignore=["u", "d", "c"]):
    """Return invariants of ``fields`` ignoring the structures in ``ignore``. The
    invariants are simplified and filtered. Doubled up operators are removed,
    including those equivalent up to index relabellings.

    Example:
        >>> invariants(L, L, Q, db, H, ignore=["u", "c", "d"])
        [L(u0, I_0)*L(u1, I_1)*Q(u2, c0, I_2)*Eps(-I_0, -I_2)*db(u3, -c1)*H(I_3)*Eps(-I_3, -I_1),
         L(u0, I_0)*L(u1, I_1)*Eps(-I_1, -I_0)*Q(u2, c0, I_2)*db(u3, -c1)*H(I_3)*Eps(-I_3, -I_2)]

    Note:
        Currently only works for SU(2) isospin structure, since
        ``remove_relabellings`` mistakenly thinks

            (L^i L^k) (L^j eb) H^l eps(ij) eps(kl)
                            and
            (L^i L^j) (L^k eb) H^l eps(ij) eps(kl)

        are the same since they can be related by interchanging j and k.

    """
    singlets = unsimplified_invariants(*fields, ignore=ignore)
    clean_singlets = remove_relabellings(clean_operators(singlets))
    with_nice_labels = [op.fill_free_indices() for op in clean_singlets]
    return with_nice_labels
