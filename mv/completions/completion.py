#!/usr/bin/env python3

"""Functions to generate completions of operators with explicit SU(2) structure."""

from mv.tensormethod.core import IndexedField
from utils import flatten, chunks
from core import Completion, EffectiveOperator
from operators import EFF_OPERATORS
from topologies import get_topology_data
from typing import Tuple, List
import networkx as nx
from copy import deepcopy

from itertools import permutations
from sympy.utilities.iterables import multiset_partitions

from functools import lru_cache


@lru_cache(maxsize=None)
def replace(data, to_replace, replace_with, found=False) -> Tuple[tuple, bool]:
    """Replace first occurance of ``to_replace`` with ``replace_with`` in
    ``data``.

    Example:
        >>> replace((("F", 18), ("S", 45), ...), "F", L('u1 i1'))
        ((L(u1, i1), 18), ("S", 45), ...), True

    """
    if found:
        return data, found

    if isinstance(data, tuple):
        new_data = []
        for datai in data:
            new_datai, found = replace(datai, to_replace, replace_with, found)
            new_data.append(new_datai)
        return tuple(new_data), found

    if data == to_replace:
        return replace_with, True

    return data, found


def replace_fields(fields, partition):
    """Takes the fields and puts them in place of the strings in the partition
    template.

        >>> replace_fields([H('i0_'), H('i1_'), L('u0_ i2_'), L('u1_ i3_')], (('F', 18), ('S', 162), (('F', 6), ('S', 54))))
        ((L(u0_, i2_), 18), (H(i0_), 162), ((L(u1_, i3_), 6), (H(i1_), 54)))

    """
    for field in fields:
        char = "S" if field.is_scalar else "F"
        partition, _ = replace(data=partition, to_replace=char, replace_with=field)

    return partition


def quick_remove_equivalent_partitions(partitions):
    # TODO Make this nicer
    return list(set(partitions))


def distribute_fields(fields, partition):
    """Takes the fields and puts them in place of the strings in the partition
    template in every possible way.

        >>> distribute_fields([H('i0_'), H('i1_'), L('u0_ i2_'), L('u1_ i3_')], (('F', 18), ('S', 162), (('F', 6), ('S', 54))))
        [((L(u0_, i2_), 18), (H(i0_), 162), ...), ((L(u1_, i3_), 18), (H(i0_), 162), ...), ...]

    Returns lots of double ups.

    """
    perms = permutations(fields)
    partitions = [replace_fields(fields, partition) for fields in perms]
    return quick_remove_equivalent_partitions(partitions)


def node_dictionary(partition):
    """Returns a dictionary mapping node to indexed field.

    Example:
        >>> node_dictionary((((Q(u364_, c210_, i369_), 6), (L(u362_, i367_), 18)),
                            ((L(u361_, i366_), 54), (Q(u363_, c209_, i368_), 162)),
                            ((db(u368_, -c214_), 486), (db(u366_, -c212_), 1458))))
        {6: Q(u364_, c210_, i369_), 18: L(u362_, i367_), ...}

    """
    flat_data = list(flatten(partition))
    tuples = chunks(flat_data, 2)
    reversed_data = map(reversed, tuples)
    return {k: {"external": v} for k, v in reversed_data}


def set_external_fields(partition, graph):
    """Add indexed fields as node attributes on graph through side effect."""
    g = deepcopy(graph)
    attrs = node_dictionary(partition)
    nx.set_node_attributes(g, attrs)
    return g


def partitions(op: EffectiveOperator) -> List[dict]:
    """Returns a list of operator partitions, epsilons and graphs of the form:

    {"fields": ((L(u0, I_0), 18), ...)
    "epsilons": (...),
    "graph": ...}

    from the partitions of the fields in the operator. This is all of the
    information required to find the completion.

    """
    topology_data_list = get_topology_data(**op.topology_type)

    out = []
    counter = 1
    for topology_data in topology_data_list:
        print(f"Furnishing topology {counter}...")

        fields = op.indexed_fields
        epsilons = op.operator.epsilons
        perms = distribute_fields(fields, topology_data["partition"])

        for perm in perms:
            g = topology_data["graph"]
            g = set_external_fields(perm, g)

            data = {"partition": perm, "epsilons": epsilons, "graph": g}
            out.append(data)

        counter += 1

    return out


def are_equivalent_partitions(a, b):
    ga = a["graph"]
    gb = b["graph"]
    node_matcher = nx.algorithms.isomorphism.categorical_node_match(
        ["external"], ["external", -1]
    )
    graph_matcher = nx.algorithms.isomorphism.GraphMatcher(
        ga, gb, node_match=node_matcher
    )
    return graph_matcher.is_isomorphic()


def remove_isomorphic(partitions: List[dict]) -> List[dict]:
    """Same algorithm as removeIsomorphic in ``wolfram/`` directory. Remove
    isomorphic graphs to reduce double-ups of completions.

    """
    parts_copy = deepcopy(partitions)

    i = 0
    while i < len(parts_copy) - 1:
        j = i + 1
        while j <= len(parts_copy) - 1:
            if are_equivalent_partitions(parts_copy[i], parts_copy[j]):
                parts_copy.pop(j)
            else:
                j += 1
        i += 1

    return parts_copy


def cons_completion(partition, epsilons, graph):
    pass
