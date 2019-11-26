#!/usr/bin/env python3

"""Functions to generate completions of operators with explicit SU(2) structure."""

from core import Completion, EffectiveOperator
from operators import EFF_OPERATORS
from topologies import get_topology_data

from itertools import permutations


def replace_field(field, char, partition, out=tuple()):
    """Replace field in partition template.

        >>> replace_field(L(u1, i1), (("F", 18), ("S", 45), ...))
        ((L(u1, i1), 18), ("S", 45), ...)

    """
    is_leaf = lambda x: not isinstance(x[0], tuple)

    for expr in partition:
        if is_leaf(expr) and expr[0] == char:
            if field is not None:
                out += ((field, expr[1]),)
                field = None
            else:
                out += (expr,)

        if is_leaf(expr) and expr[0] != char:
            out += (expr,)

        if not is_leaf(expr):
            out += (replace_field(field, char, expr),)

    return out


def replace_fields(fields, partition):
    """Takes the fields and puts them in place of the strings in the partition
    template.

        >>> replace_fields([H(i0_), H(i1_), L(u0_, i2_), L(u1_, i3_)], (('F', 18), ('S', 162), (('F', 6), ('S', 54))))
        ((L(u0_, i2_), 18), (H(i0_), 162), ...)

    """
    for field in fields:
        char = "S" if field.is_scalar else "F"
        partition = replace_field(field, char, partition)

    return partition


def remove_equivalent_partitions(partitions):
    # TODO Make this nicer
    return list(set(partitions))


def distribute_fields(fields, partition):
    """Takes the fields and puts them in place of the strings in the partition
    template in every possible way.

        >>> distribute_fields([H(i0_), H(i1_), L(u0_, i2_), L(u1_, i3_)], (('F', 18), ('S', 162), (('F', 6), ('S', 54))))
        [((L(u0_, i2_), 18), (H(i0_), 162), ...), ((L(u1_, i3_), 18), (H(i0_), 162), ...), ...]

    Returns lots of double ups.

    """
    perms = permutations(fields)
    partitions = [replace_fields(fields, partition) for fields in perms]
    return remove_equivalent_partitions(partitions)


def operator_partition(op: EffectiveOperator, partition) -> dict:
    """Returns an operator partition of the form:

    {"fields": ((L(u0, I_0), 18), ...)
    "epsilons": (...)}

    from the partition template.

    """
    fields = op.indexed_fields
    epsilons = op.operator.structures

    perms = distribute_fields(fields, partition)


def get_completions(op: EffectiveOperator):
    topology_data = get_topology_data(**op.topology_type)
    completions = [Completion(**topology) for topology in topolopy_data]
    field_perms = permutations(op.indexed_fields)
    pass
