#!/usr/bin/env python3

"""Functions to generate completions of operators with explicit SU(2) structure."""

from mv.tensormethod.core import (
    Index,
    IndexedField,
    eps,
    delta,
    is_invariant_symbol,
    Operator,
    get_dynkin,
)
from mv.tensormethod.contract import (
    lorentz_singlets,
    colour_singlets,
    invariants,
    contract_su2,
)
from utils import flatten, chunks
from core import (
    Completion,
    FailedCompletion,
    EffectiveOperator,
    cons_completion_field,
    FieldType,
    VectorLikeDiracFermion,
    MajoranaFermion,
    ComplexScalar,
)
from operators import EFF_OPERATORS
from topologies import get_topology_data, Leaf
from typing import Tuple, List, Dict, Union
import networkx as nx
from copy import deepcopy

from itertools import permutations, groupby
from sympy.utilities.iterables import multiset_partitions
from sympy.tensor.tensor import Tensor

from functools import lru_cache, reduce


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


def _fix_su2_indices(free_indices: List[Index]) -> Tuple[list, list]:
    """UNUSED NOW

    Raises lowered SU2 indices in free_indces and returns new epsilons.

    """
    new_epsilons, new_free_indices = [], []
    for idx in free_indices:
        if idx.is_su2 and not idx.is_up:
            type_ = Index.get_index_labels()[idx.index_type]
            new_index = Index.fresh(type_)
            new_free_indices.append(new_index)
            new_epsilons.append(eps(str(-new_index) + f" {idx}"))
        else:
            new_free_indices.append(idx)

    return new_free_indices, new_epsilons


def is_contracted_epsilon(eps, indices):
    i, j, *k = eps.indices

    # deal with su2 epsilon first
    if not k:
        if -i in indices and -j in indices:
            indices.remove(-i)
            indices.remove(-j)
            return True
        return False

    # su3 epsilon has three indices
    else:
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
) -> Tuple[FieldType, List[Tensor], List[Operator]]:
    """Takes two or three indexed fields and the epsilons and returns a new indexed field
    transforming in the same way as x \otimes y.

    Returns a possibly empty tuple.

    Example:
        >>> s = contract((H('i0'), H('i1')), 'S', epsilons=[])
        >>> s
        S(i0, i1)
        >>> s.y
        1

    """
    if len(fields) != 2 and len(fields) != 3:
        raise Exception("Too many fields passed to contract.")

    # product of fields
    prod_fields = reduce(lambda a, b: a * b, fields)
    free_indices = prod_fields.get_free_indices()

    # remove contracted isospin epsilons and remove indices from free_indices by
    # side effect
    contracted_epsilons, spectator_epsilons = [], []
    # contracted epsilons are those all of whose indices are contracted on the
    # set of fields passed in. Spectator epsilons are all of the others
    for e in epsilons:
        if is_contracted_epsilon(e, free_indices):
            contracted_epsilons.append(e)
        else:
            spectator_epsilons.append(e)

    undotted, dotted, colour, isospin, _ = Index.indices_by_type(free_indices).values()

    # sort out lorentz epsilons from contraction
    scalar = len(undotted) == 0 and len(dotted) == 0
    two_undotted = len(undotted) == 2 and len(dotted) == 0
    two_dotted = len(dotted) == 2 and len(undotted) == 0
    one_undotted = len(undotted) == 1 and len(dotted) == 0
    one_dotted = len(dotted) == 1 and len(undotted) == 0

    if not (scalar or two_undotted or two_dotted or one_undotted or one_dotted):
        return None

    # construct new (lorentz) epsilons for scalar contractions (currently not
    # allowing vector contractions)
    lorentz_contraction_epsilons = []
    if two_undotted:
        eps_indices = " ".join(str(-i) for i in undotted)
        lorentz_eps = eps(eps_indices)
        lorentz_contraction_epsilons.append(lorentz_eps)
        # remove newly contracted indices
        for i in undotted:
            free_indices.remove(i)
    elif two_dotted:
        eps_indices = " ".join(str(-i) for i in dotted)
        lorentz_eps = eps(eps_indices)
        lorentz_contraction_epsilons.append(lorentz_eps)
        # remove newly contracted indices
        for i in dotted:
            free_indices.remove(i)

    # construct exotic field
    # deal with charges
    exotic_charges = {}
    pairs = map(lambda x: x.charges.items(), fields)
    for n_pairs in zip(*pairs):
        # fish out key
        (k, _), *rst = n_pairs
        # sum over values
        exotic_charges[k] = sum(map(lambda x: x[1], n_pairs))

    indices_by_type = Index.indices_by_type(free_indices).values()
    exotic_undotted, exotic_dotted, exotic_colour, exotic_isospin, _ = map(
        sorted, indices_by_type
    )

    exotic_indices = " ".join(
        str(i)
        for i in [*exotic_undotted, *exotic_dotted, *exotic_colour, *exotic_isospin]
    )

    # establish fermion or boson for symbols
    lorentz_dynkin = get_dynkin(exotic_indices)[:2]
    if lorentz_dynkin == "10" or lorentz_dynkin == "01":
        symbols_to_use = symbols["fermion"]
    else:
        symbols_to_use = symbols["boson"]

    fs = sorted([f.field for f in fields], key=lambda x: x.label_with_dagger)
    fs = tuple(fs)
    if fs in field_dict.keys():
        symbol = field_dict[fs]
    else:
        symbol = symbols_to_use.pop(0)
        field_dict[fs] = symbol

    # mutate available symbols
    exotic_field = IndexedField(symbol, exotic_indices, charges=exotic_charges)

    # construct MajoranaFermion, VectorLikeDiracFermion, ...
    exotic_field = cons_completion_field(exotic_field)
    partner = exotic_field
    if isinstance(exotic_field, VectorLikeDiracFermion):
        partner = exotic_field.dirac_partner()
    elif isinstance(exotic_field, ComplexScalar):
        partner = exotic_field.conj
    elif isinstance(exotic_field, MajoranaFermion):
        partner = exotic_field.majorana_partner()

    # need additional su2 epsilons to fix su2 indices (since not working with
    # lowered indices at all)
    partner, fix_su2_epsilons = partner.lower_su2()

    # prod_epsilons = new_epsilons + contracted_epsilons
    # prod = reduce(lambda a, b: a * b, prod_epsilons, prod_fields)

    # construct term
    new_epsilons = [
        *fix_su2_epsilons,
        *contracted_epsilons,
        *lorentz_contraction_epsilons,
    ]
    term = reduce(lambda a, b: a * b, new_epsilons, prod_fields)
    term *= partner
    assert not term.free_indices

    return exotic_field, term, spectator_epsilons, lorentz_contraction_epsilons


def get_connecting_edge(graph: nx.Graph, nodes: List[int]) -> Tuple[int, int]:
    """Returns an edge that connects to nodes in ``nodes``.

    Taking ``graph`` to be:

            4
            |
            |
            3
           / \
          /   \
         1     2

    Example:
        >>> get_connecting_edge(graph, (1, 2))
        (3, 4)

    """
    neighbours = {}
    for node in nodes:
        neighbours[node] = set(graph.neighbors(node))

    fst, *rst = list(neighbours.values())
    intersection = fst.intersection(*rst)
    assert len(intersection) == 1
    connecting_node = list(intersection)[0]

    other_nodes = list(graph.neighbors(connecting_node))
    for node in nodes:
        other_nodes.remove(node)

    assert len(other_nodes) == 1
    return (connecting_node, other_nodes[0])


def build_term(
    leaves: Tuple[Tuple[IndexedField, int]],
    symbols: List[str],
    epsilons: list,
    lorentz_epsilons: list,
    terms: list,
    edge_dict: Dict[FieldType, Tuple[int, int]],
    field_dict: Dict[tuple, str],
    graph: nx.Graph,
) -> Tuple[IndexedField, Tuple[int, int]]:
    """Takes two leaf-tuples of indexed fields and vertex labels, constructs
    lagrangian term and appends it to ``terms`` and returns new leaf-tuple like
    (exotic indexed field, corresponding edge).

    Updates terms, edge_dict and epsilons by side effect.

    Example:
        >>> terms, edges, epsilons = [], {}, [metric(-i366_, -i368_), metric(-i369_, -i367_)]
        >>> build_term(((L(u362_, i367_), 6), (L(u361_, i366_), 18)),
                       epsilons=epsilons,
                       terms=terms,
                       edge_dict=edges,
                       graph=G)

    """
    # update epsilons by side effect
    fields, nodes = [], []
    for leaf in leaves:
        if leaf == (None, None):
            return Leaf(None, None)

        field, node = leaf
        fields.append(field)
        nodes.append(node)

    # if only one field, SM field at last vertex
    if len(fields) == 1:
        return Leaf(field, node)

    # update epsilons
    try_contract = contract(fields, symbols, epsilons, field_dict)

    if try_contract is None:
        return Leaf(None, None)

    exotic_field, term, epsilons, new_epsilons = try_contract
    lorentz_epsilons += new_epsilons

    exotic_edge = get_connecting_edge(graph, nodes)

    # update edge_dict
    # keep track of which fields go where in the diagram
    edge_dict[exotic_field] = exotic_edge

    # update terms
    terms.append(term)
    return Leaf(exotic_field, exotic_edge[0])


def contains_only_leaves(xs):
    if not isinstance(xs, tuple):
        return False

    for x in xs:
        if not isinstance(x, Leaf):
            return False

    return True


def reduced_row(row, func):
    if isinstance(row, Leaf):
        return row

    # row is a tuple
    if contains_only_leaves(row):
        return func(row)

    return func(tuple(map(lambda a: reduced_row(a, func), row)))


def cons_completion(
    partition, epsilons, graph, filter_function=None
) -> Tuple[str, tuple]:
    """Work in progress...

    ``filter_function`` is a dyadic function that takes an IndexedField and an
    interaction operator and returns a bool. Implement no vectors in completions
    like

        >>> filter_func = lambda f, int: not f.is_vector

    """
    lorentz_epsilons, terms, edge_dict, field_dict = [], [], {}, {}
    symbols = {"fermion": ["ψ", "ξ", "χ", "ζ", "f"], "boson": ["φ", "η", "ω", "ρ", "S"]}

    func = lambda leaves: build_term(
        leaves=leaves,
        symbols=symbols,
        epsilons=epsilons,
        lorentz_epsilons=lorentz_epsilons,
        terms=terms,
        edge_dict=edge_dict,
        field_dict=field_dict,
        graph=graph,
    )

    reduced_partition = [reduced_row(row, func) for row in partition]

    # Add last interaction to term
    prod = None
    for i in reduced_partition:
        f = i.field
        if f is None:
            return None

        if prod is None:
            prod = f
        else:
            prod *= f

        # if isinstance(f, FieldType):
        #     prod *= f.args[0](*f.args[1])
        # else:
        #     prod *= f

    free_indices = prod.get_free_indices()
    undotted, dotted, colour, isospin, _ = Index.indices_by_type(free_indices).values()

    if len(undotted) == 2 and len(dotted) == 0:
        eps_indices = " ".join(str(-i) for i in undotted)
        lorentz_epsilons.append(eps(eps_indices))
        prod *= eps(eps_indices)
    elif len(undotted) == 0 and len(dotted) == 2:
        eps_indices = " ".join(str(-i) for i in dotted)
        lorentz_epsilons.append(eps(eps_indices))
        prod *= eps(eps_indices)
    elif not (len(undotted) == 0 and len(dotted) == 0):
        return None

    free_indices = prod.free_indices
    contracted_epsilons, spectator_epsilons = [], []
    for e in epsilons:
        if is_contracted_epsilon(e, free_indices):
            contracted_epsilons.append(e)
        else:
            spectator_epsilons.append(e)

    prod = reduce(lambda a, b: a * b, contracted_epsilons, prod)
    assert not prod.free_indices

    terms.append(prod)

    return terms, edge_dict, field_dict, lorentz_epsilons


def partition_completion(partition) -> Union[Completion, FailedCompletion]:
    part = partition["partition"]
    epsilons = partition["epsilons"]
    graph = partition["graph"]
    op = partition["operator"]

    # failed is a string with the reason the completion failed
    try:
        args = cons_completion(partition=part, epsilons=epsilons, graph=graph)
    except:
        breakpoint()
    if args is not None:
        terms, edge_dict, field_dict, lorentz_epsilons = args
    else:
        return FailedCompletion("Bad lorentz irrep")

    prod = reduce(lambda a, b: a * b, lorentz_epsilons, 1)
    explicit_op = op.operator * prod
    exotics = set(f for f in edge_dict.keys())
    eff_operator = EffectiveOperator(op.name, Operator(explicit_op))

    return Completion(
        operator=eff_operator, partition=part, graph=graph, exotics=exotics, terms=terms
    )


def operator_completions(operator) -> List[Completion]:
    parts = partitions(operator)
    completions = [partition_completion(p) for p in parts]
    good_completions = [c for c in completions if not isinstance(c, FailedCompletion)]
    return good_completions


def lnv_completions(operator_name: str) -> List[Completion]:
    op = EFF_OPERATORS[operator_name]
    return operator_completions(op)


def collect_completions(
    completions: List[Completion]
) -> Dict["FieldContent", Completion]:
    """Return dictionary mapping field content to list of completions."""
    out = {}
    for k, g in groupby(completions, key=lambda x: x.exotic_info()):
        g = list(g)
        out[k] = g

    return out


def filter_completions(
    completions: List[Completion], sieve: List[Completion]
) -> List[Completion]:
    pass
