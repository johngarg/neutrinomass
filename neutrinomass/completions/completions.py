#!/usr/bin/env python3

"""Functions to generate completions of operators with explicit SU(2) structure."""

from neutrinomass.tensormethod.core import (
    Index,
    IndexedField,
    eps,
    delta,
    is_invariant_symbol,
    Operator,
    get_dynkin,
)
from neutrinomass.tensormethod.contract import (
    lorentz_singlets,
    colour_singlets,
    invariants,
    term_lorentz_singlet,
)

from neutrinomass.tensormethod.utils import safe_nocoeff
from neutrinomass.completions.utils import (
    flatten,
    chunks,
    factors,
    remove_equivalent,
    multiple_replace,
)
from neutrinomass.completions.core import (
    Completion,
    FailedCompletion,
    EffectiveOperator,
    cons_completion_field,
    FieldType,
    VectorLikeDiracFermion,
    MajoranaFermion,
    ComplexScalar,
)
from neutrinomass.completions.topologies import get_topology_data, Leaf
from neutrinomass.utils import pmatch

from typing import Tuple, List, Dict, Union
import networkx as nx
from copy import deepcopy

from itertools import permutations, groupby
from sympy.tensor.tensor import Tensor
from sympy import prime

from functools import lru_cache, reduce
import re


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

        f = lambda x: Leaf(*x) if isinstance(data, Leaf) else tuple(x)
        return f(new_data), found

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
        char = "S" if field.is_boson else "F"
        partition, _ = replace(data=partition, to_replace=char, replace_with=field)

    return partition


def quick_remove_equivalent_partitions(partitions):
    """Just remove double ups. (For now.)

    """
    return list(set(partitions))


def distribute_fields(fields, partition):
    """Takes the fields and puts them in place of the strings in the partition
    template in every possible way.

        >>> distribute_fields([H('i0_'), H('i1_'), L('u0_ i2_'), L('u1_ i3_')], (('F', 18), ('S', 162), (('F', 6), ('S', 54))))
        [((L(u0_, i2_), 18), (H(i0_), 162), ...), ((L(u1_, i3_), 18), (H(i0_), 162), ...), ...]

    Returns lots of double ups.

    """
    perms = permutations(fields)
    parts = [replace_fields(fields, partition) for fields in perms]
    return quick_remove_equivalent_partitions(parts)


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
    return {k: {"particle": v} for k, v in reversed_data}


def set_external_fields(partition, graph):
    """Add indexed fields as edge attributes on graph through side effect."""
    g = deepcopy(graph)
    node_attrs = node_dictionary(partition)

    edge_attrs = {}
    for edge in graph.edges:
        for n, field_dict in node_attrs.items():
            if n in edge:
                edge_attrs[edge] = field_dict

    nx.set_edge_attributes(g, edge_attrs)
    return g


def partitions(operator: EffectiveOperator) -> List[dict]:
    """Returns a list of operator partitions, epsilons and graphs of the form:

    {"fields": ((L(u0, I_0), 18), ...)
    "epsilons": (...),
    "graph": ...}

    from the partitions of the fields in the operator. This is all of the
    information required to find the completion.

    """
    topology_data_list = get_topology_data(**operator.topology_type)

    colour_ops = colour_singlets([operator.operator], overcomplete=True)
    colour_ops = [EffectiveOperator(operator.name, op) for op in colour_ops]

    out = []
    counter = 1
    for topology_data in topology_data_list:
        # print(f"Furnishing topology {counter}...")
        for op in colour_ops:
            fields = op.indexed_fields
            epsilons = op.operator.epsilons
            perms = distribute_fields(fields, topology_data["partition"])

            for perm in perms:
                g = topology_data["graph"]
                g = set_external_fields(perm, g)

                data = {
                    "operator": op,
                    "partition": perm,
                    "epsilons": epsilons,
                    "graph": g,
                }
                out.append(data)

        counter += 1

    return out


def are_equivalent_partitions(a, b):
    ga = a["graph"]
    gb = b["graph"]
    node_matcher = nx.algorithms.isomorphism.categorical_node_match(
        ["external"], ["external", -1]  # set default value to -1 for internal nodes
    )
    graph_matcher = nx.algorithms.isomorphism.GraphMatcher(
        ga, gb, node_match=node_matcher
    )

    return graph_matcher.is_isomorphic()


def remove_isomorphic(partitions: List[dict]) -> List[dict]:
    """Same algorithm as removeIsomorphic in ``wolfram/`` directory. Remove
    isomorphic graphs to reduce double-ups of completions.

    """
    remove_equivalent(partitions, are_equivalent_partitions)
    return None


def is_contracted_epsilon(eps, indices):
    """Return True if two indices on epsilon are contracted, False otherwise."""
    i, j, *k = eps.indices

    # deal with su2 epsilon first
    if not k:
        if -i in indices and -j in indices:
            indices.remove(-i)
            indices.remove(-j)
            return True
        return False

    # su3 epsilon has three indices
    counter = 0
    to_remove, free = [], []
    for idx in eps.indices:
        if -idx in indices:
            to_remove.append(-idx)
        else:
            free.append(idx)

    # is contracted epsilon
    if len(to_remove) == 2:
        for idx in to_remove:
            indices.remove(idx)

        assert len(free) == 1
        indices.append(free[0])
        return True

    if len(to_remove) == 3:
        for idx in to_remove:
            indices.remove(idx)

        return True

    return False


def contract(
    fields: Tuple[IndexedField],
    symbols: Dict[str, List[str]],
    epsilons: list,
    field_dict: Dict[tuple, str],
) -> Tuple[FieldType, Operator, List[Tensor], List[Tensor]]:
    """Takes two or three indexed fields and the epsilons [epsilons and deltas of
    SU(2) and SU(3) from the operator] and returns a new indexed field
    transforming in the same way as $x \otimes y$.

    Gauge epsilons and deltas are going to be potentially used up in this
    process, while epsilons carrying Lorentz indices will be introduced
    enforcing the contractions between dotted and undotted indices in the
    generated operator.

    Returns a tuple with the field transforming like the product of `fields`,

    Example:
        >>> field, term, epsilons, new_epsilons = contract((H('i0'), H('i1')), {"fermion": [], "boson": ["S"]}, [], {})
        >>> field
        S(i0, i1)
        >>> field.y
        1

    """
    pass
