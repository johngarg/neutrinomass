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
    contract_su2,
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
from copy import copy, deepcopy

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


# The approach to finding the completions is the following: contract off fields
# and find corresponding exotic and term. Replace the fields by the exotic and
# keep track of the available epsilons and the terms by mutation. The pipeline is
#
# contract: returns exotic field, new gauge epsilons (fewer) and new lorentz
# epsilons (more)
#
# replace_and_mutate: returns a Leaf structure that enters the partition in
# place of the contracted fields, mutates terms, edge_dict of graph,
# gauge_epsilons and lorentz_epsilons
#
# reduce_partition: applies replace_and_mutate to a partition until last vertex.


def all_scalars(fields):
    boolean = True
    for f in fields:
        boolean = boolean and f.is_boson

    return boolean


def all_fermions(fields):
    boolean = True
    for f in fields:
        boolean = boolean and f.is_fermion

    return boolean


def drop_scalar(fields):
    scalars, fermions = [], []
    for f in fields:
        if f.is_boson:
            scalars.append(f)
        elif f.is_fermion:
            fermions.append(f)
    assert len(scalars) == 1
    return fermions


def get_lorentz_epsilons(fields: Tuple[IndexedField]) -> Tuple[bool, List[Tensor]]:
    """Takes a list of two or three fields (possibly with derivatives) and returns
    the lorentz epsilons that contract the fields to as low a Lorentz irrep as
    possible as well as a boolean indicating whether the contraction is allowed.

    """

    deriv_structure = [f.derivs for f in fields]
    n_derivs = sum(deriv_structure)

    if n_derivs > 2:
        raise Exception(
            f"Not currently supporting {n_derivs} derivatives in an operator."
        )

    if not n_derivs and len(fields) == 4:
        return True, []

    if not n_derivs and len(fields) == 3:
        if all_scalars(fields):
            return True, []

        elif all_fermions(fields):
            return False, []

        return get_lorentz_epsilons(drop_scalar(fields))

    if n_derivs == 2 and len(fields) == 3:
        fields = sorted(fields, key=lambda f: -f.derivs)

    prod = reduce(lambda x, y: x * y, fields)
    undotted, dotted, _, _, _, = prod.indices_by_type.values()

    # Reject vector contraction
    if len(undotted) == 1 and len(dotted) == 1:
        return False, []

    epsilons = []
    for indices in [undotted, dotted]:
        # skip single indices (fermion, scalar) contraction
        if len(indices) == 1:
            continue

        # pair up all even indices; if odd, leave last index
        if len(indices) % 2 != 0:
            indices.pop(-1)

        for i, j in chunks(indices, 2):
            epsilons.append(eps(f"-{i} -{j}"))

    return True, epsilons


def is_contracted_epsilon(eps, indices):
    """Return True if two indices on epsilon are contracted, False otherwise."""
    i, j, *k = eps.indices

    # deal with su2 epsilon first
    if not k:
        if -i in indices and -j in indices:
            return True
        return False

    # su3 epsilon has three indices
    to_remove, free = [], []
    for idx in eps.indices:
        if -idx in indices:
            to_remove.append(-idx)
        else:
            free.append(idx)

    # is contracted epsilon
    if len(to_remove) == 2:
        assert len(free) == 1
        indices.append(free[0])
        return True

    if len(to_remove) == 3:
        return True

    return False


def separate_gauge_epsilons(
    fields: List[IndexedField], epsilons: List[Tensor]
) -> Tuple[List[Tensor], List[Tensor]]:

    prod_fields = reduce(lambda x, y: x * y, fields)
    free_indices = prod_fields.get_free_indices()
    contracted_epsilons, spectator_epsilons = [], []

    # contracted epsilons are those all of whose indices are contracted on the
    # set of fields passed in. Spectator epsilons are all of the others
    for epsilon in epsilons:
        if is_contracted_epsilon(epsilon, free_indices):
            contracted_epsilons.append(epsilon)
        else:
            spectator_epsilons.append(epsilon)

    return spectator_epsilons, contracted_epsilons


def check_charges(operator, ignore=[]):
    fields = [f for f in operator.tensors if isinstance(f, IndexedField)]
    for q in fields[0].charges:
        if q in ignore:
            continue
        assert not sum(f.charges[q] for f in fields)


def check_singlet(operator, ignore=["3b"]):
    check_charges(operator, ignore=ignore)
    for free in operator.free_indices:
        assert free.index_type == "Generation"


def exotic_field_and_term(
    op: Operator, symbols: Dict[str, List[str]], field_dict: Dict[tuple, str]
) -> Tuple[IndexedField, IndexedField, Union[Operator, str]]:
    """Returns exotic field, partner (that couples in Lagrangian) and Lagrangian
    term. Mutates field_dict with exotic field's symbol.

    The last item returned may be a string explaining why the contraction
    failed.

    """

    fields = [f for f in op.tensors if isinstance(f, IndexedField)]
    exotic_charges = {}
    pairs = list(map(lambda x: x.charges.items(), fields))

    # ensure no fields have missing or extra charges
    fst, *rst = pairs
    for pair in rst:
        assert len(pair) == len(fst)

    for n_pairs in zip(*pairs):
        # fish out key
        (k, _), *rst = n_pairs
        # sum over values
        exotic_charges[k] = sum(map(lambda x: x[1], n_pairs))

    indices_by_type = op.indices_by_type.values()
    exotic_undotted, exotic_dotted, exotic_colour, exotic_isospin, _ = map(
        sorted, indices_by_type
    )

    exotic_indices = " ".join(
        str(i)
        for i in [*exotic_undotted, *exotic_dotted, *exotic_colour, *exotic_isospin]
    )

    # establish fermion or boson for symbols
    lorentz_dynkin = get_dynkin(exotic_indices)[:2]
    if lorentz_dynkin in ("10", "01"):
        symbols_to_use = symbols["fermion"]
    else:
        symbols_to_use = symbols["boson"]

    fs = sorted([f.field for f in fields], key=lambda x: x.label_with_dagger)
    fs = tuple(fs)
    if fs in field_dict.keys():
        symbol = field_dict[fs]
    else:
        symbol = symbols_to_use.pop(0)
        # field_dict mutated here!
        field_dict[fs] = symbol

    # for Dirac and Majorana fermions, always keep plain symbol left handed
    to_conj = False
    if exotic_indices and exotic_indices[0] == "d":
        exotic_indices = Index.conj_index_string(exotic_indices)
        to_conj = True

    exotic_indexed_field = IndexedField(
        label=symbol, indices=exotic_indices, charges=exotic_charges
    )

    # construct MajoranaFermion, VectorLikeDiracFermion, ...
    # `partner` is in the Lagrangian, `exotic_field` transforms like contracted pair
    exotic_field = cons_completion_field(exotic_indexed_field)
    exotic_field = exotic_field.conj_indices if to_conj else exotic_field

    partner = exotic_field
    if isinstance(exotic_field, VectorLikeDiracFermion):
        partner = exotic_field.dirac_partner()
    elif isinstance(exotic_field, ComplexScalar):
        partner = exotic_field.conj
    elif isinstance(exotic_field, MajoranaFermion):
        partner = exotic_field.majorana_partner()

    # Need additional su2 epsilons to fix su2 indices (since not working with
    # lowered indices at all). Won't need to do this if removing a derivative in
    # the process
    partner, fix_su2_epsilons = partner.lower_su2()
    term = reduce(lambda x, y: x * y, fix_su2_epsilons, op * partner)

    # construct term and check to see if vanishes
    if term.canon_bp() == 0:
        return exotic_field, partner, f"Vanishing coupling at {term}"

    # need to construct term again because sympy is annoying
    term = reduce(lambda x, y: x * y, fix_su2_epsilons, op * partner)

    check_singlet(term)

    return exotic_field, partner, term


def process_derivative_term(op: Operator) -> Union[Operator, str]:
    deriv_structure = [f.derivs for f in op.fields]
    n_derivs = sum(deriv_structure)

    if n_derivs == 0:
        return op

    # Remove derivatives and lorentz epsilons and call contract_su2 on Lorentz
    # structure.
    #
    # There are a number of cases here
    # 1. (DH)(DH)SS
    # 2. (DH)(DH)S
    # 3. (DH)ψF -> in all but this case, affected fields are clear
    # 4. S(Dψ)F
    # 5. (Dψ)(Dψ)S
    scalars, fermions, exotic_fermions, epsilons = [], [], [], []
    for t in op.tensors:
        if isinstance(t, IndexedField) and t.derivs > 0 and t.is_boson:
            no_deriv_field = t.strip_derivs_with_indices()
            scalars.append(no_deriv_field)
        if isinstance(t, IndexedField) and t.derivs > 0 and t.is_fermion:
            no_deriv_field = t.strip_derivs_with_indices()
            fermions.append(no_deriv_field)
        elif isinstance(t, VectorLikeDiracFermion):
            fixed_field = t.dirac_partner().conj
            exotic_fermions.append(fixed_field)
        elif isinstance(t, MajoranaFermion):
            fixed_field = t.conj
            exotic_fermions.append(fixed_field)
        elif isinstance(t, IndexedField) and t.is_fermion:
            fermions.append(t)
        elif isinstance(t, IndexedField) and t.is_scalar:
            scalars.append(t)
        elif not isinstance(t, IndexedField):
            # is epsilon, keep gauge ones, not lorentz
            if t.indices[0].index_type in ("Undotted", "Dotted"):
                continue
            else:
                epsilons.append(t)

    # cases 1 and 2
    if len(scalars) > 2:
        return reduce(lambda x, y: x * y, scalars + epsilons)
    # case 5
    if len(fermions) == 2:
        left, right = fermions
    # cases 3 and 4
    if len(exotic_fermions) == 1:
        assert len(fermions) == 1
        left, right = exotic_fermions[0], fermions[0]
    if len(exotic_fermions) == 2:
        left, right = exotic_fermions

    # include scalars and epsilons in su2 contraction
    right = reduce(lambda x, y: x * y, scalars + epsilons, right)
    lu, ld, _, _, _ = left.indices_by_type.values()
    ru, rd, _, _, _ = right.indices_by_type.values()

    # if the indices are equal after taking the conj, then there will be an
    # error. In this case, you can just lower one of them
    if lu == ru and ld == rd:
        partner, fix_su2_epsilons = left.lower_su2(skip=["Isospin"])
        assert len(fix_su2_epsilons) == 1
        return right * partner * fix_su2_epsilons[0]
    if lu:
        if not (len(lu) == 1 and len(ru) == 1):
            return "Not allowed contraction"
        index_str = " ".join(str(-i) for i in lu + ru)
    else:
        if not (len(ld) == 1 and len(rd) == 1):
            return "Not allowed contraction"
        index_str = " ".join(str(-i) for i in ld + rd)

    return right * left * eps(index_str)


def contract(
    fields: Tuple[IndexedField],
    symbols: Dict[str, List[str]],
    gauge_epsilons: list,
    field_dict: Dict[tuple, str],
) -> Union[Tuple[FieldType, Operator, List[Tensor], List[Tensor]], str]:
    """Takes two or three indexed fields and the epsilons [epsilons and deltas of
    SU(2) and SU(3) from the operator] and returns a new indexed field
    transforming in the same way as $x \otimes y$.

    Gauge epsilons (and deltas) are going to be potentially used up in this
    process, while epsilons carrying Lorentz indices will be introduced
    enforcing the contractions between dotted and undotted indices in the
    generated operator.

    Returns a tuple with the field transforming like the product of `fields`,
    the term, and the new gauge and lorentz epsilons.

    If the contraction fails, returns a string with the reason it failed.

    Example:
        >>> field, term, gauge_epsilons, lorentz_epsilons = contract((H('i0'), H('i1')), [], {"fermion": [], "boson": ["S"]}, {})
        >>> field
        S(i0, i1)
        >>> field.y
        1

    """
    if len(fields) != 2 and len(fields) != 3:
        raise Exception("Too many fields passed to contract.")

    allowed_contraction, lorentz_epsilons = get_lorentz_epsilons(fields)
    if not allowed_contraction:
        # Bad lorentz contraction
        return "Bad Lorentz contraction."

    # some gauge epsilons will be removed in the contraction, the others will
    # just watch
    spectator_gauge_eps, eps_to_remove = separate_gauge_epsilons(fields, gauge_epsilons)
    # contracted_fields is the fields (with the derivatives still present) with
    # lorentz indices contracted
    contracted_fields = reduce(
        lambda x, y: x * y, (*fields, *lorentz_epsilons, *eps_to_remove)
    )

    exotic, partner, maybe_term = exotic_field_and_term(
        contracted_fields, symbols, field_dict
    )

    if isinstance(maybe_term, str):
        # return the reason
        return maybe_term

    no_deriv_maybe_term = process_derivative_term(maybe_term)
    if isinstance(no_deriv_maybe_term, str):
        return no_deriv_maybe_term

    if no_deriv_maybe_term.canon_bp() == 0:
        return f"Vanishing coupling at {maybe_term} after derivative processing."

    return exotic, no_deriv_maybe_term, spectator_gauge_eps, lorentz_epsilons


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


def replace_and_mutate(
    leaves: Tuple[Tuple[IndexedField, int]],
    symbols: List[str],
    gauge_epsilons: list,
    lorentz_epsilons: list,
    terms: list,
    edge_dict: Dict[FieldType, Tuple[int, int]],
    field_dict: Dict[tuple, str],
    graph: nx.Graph,
) -> Tuple[IndexedField, Tuple[int, int]]:
    """Returns a Leaf structure that enters the partition in place of the contracted
    fields. Mutates major state of completion: terms, edge_dict of graph,
    gauge_epsilons and lorentz_epsilons. Mutation of the field_dict happens in
    `exotic_field_and_term` through `contract`.

    For a failed completion, keep reason in first element of leaf-tuple.

    TODO Fix example
    Example:
        >>> terms, edges, epsilons =
        >>> replace_and_mutate(((L(u362_, i367_), 6), (L(u361_, i366_), 18)),
                               symbols=symbols,
                               gauge_epsilons=gauge_epsilons,
                               lorentz_epsilons=lorentz_epsilons,
                               terms=terms,
                               edge_dict=edges,
                               field_dict=field_dict,
                               graph=G)

    """
    fields, nodes = [], []
    for leaf in leaves:
        if leaf[1] == None:
            return leaf

        field, node = leaf
        fields.append(field)
        nodes.append(node)

    # if only one field, SM field at last vertex
    if len(fields) == 1:
        return Leaf(field, node)

    # field_dict is updated in this call
    maybe_contract = contract(fields, symbols, gauge_epsilons, field_dict)

    # For a failed completion, keep reason in first element of leaf-tuple.
    if isinstance(maybe_contract, str):
        return Leaf(maybe_contract, None)

    # mutate gauge_epsilons immediately
    exotic_field, term, gauge_epsilons, new_lorentz_epsilons = maybe_contract

    # mutate lorentz_epsilons, terms
    lorentz_epsilons += new_lorentz_epsilons
    terms.append(term)

    # update edge_dict
    exotic_edge = get_connecting_edge(graph, nodes)
    edge_dict[exotic_field] = exotic_edge

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


def construct_completion(partition, gauge_epsilons, graph):
    lorentz_epsilons, terms, edge_dict, field_dict = [], [], {}, {}
    symbols = {"fermion": ["ψ", "ξ", "χ", "ζ", "f"], "boson": ["φ", "η", "ω", "ρ", "S"]}

    func = lambda leaves: replace_and_mutate(
        leaves=leaves,
        symbols=symbols,
        gauge_epsilons=gauge_epsilons,
        lorentz_epsilons=lorentz_epsilons,
        terms=terms,
        edge_dict=edge_dict,
        field_dict=field_dict,
        graph=graph,
    )

    reduced_partition = [reduced_row(row, func) for row in partition]

    # construct final interaction term and add to terms
    prod = None
    for i in reduced_partition:
        f = i.field
        if isinstance(f, str):
            return f

        if prod is None:
            prod = f
        else:
            prod *= f

    fields = [f for f in prod.tensors if isinstance(f, IndexedField)]

    allowed, new_lorentz_epsilons = get_lorentz_epsilons(fields)
    if not allowed:
        return "Bad Lorentz contraction."

    _, eps_to_remove = separate_gauge_epsilons(fields, gauge_epsilons)

    for e in [*new_lorentz_epsilons, *eps_to_remove]:
        prod *= e

    # mutate Lorentz epsilons with last contraction
    lorentz_epsilons += new_lorentz_epsilons

    proc_term = process_derivative_term(prod)

    if isinstance(proc_term, str):
        return proc_term

    if proc_term.canon_bp() == 0:
        return f"Vanishing coupling at {proc_term} after derivative processing."

    # make sure the term is a singlet
    check_singlet(proc_term)

    # append the processed term to terms
    terms.append(proc_term)

    return terms, edge_dict, field_dict, lorentz_epsilons


def partition_completion(partition) -> Union[Completion, FailedCompletion]:
    part = partition["partition"]
    gauge_epsilons = partition["epsilons"]
    graph = partition["graph"]
    op = partition["operator"]

    # if args is a string, then it's the reason the completion failed
    args = construct_completion(part, gauge_epsilons, graph)
    if not isinstance(args, str):
        terms, edge_dict, field_dict, lorentz_epsilons = args
    else:
        return FailedCompletion(args)

    explicit_op = reduce(lambda a, b: a * b, lorentz_epsilons, op.operator)
    exotics = set(f for f in edge_dict.keys())
    eff_operator = EffectiveOperator(op.name, explicit_op)

    new_edge_attrs = {v: {"particle": k} for k, v in edge_dict.items()}
    nx.set_edge_attributes(graph, new_edge_attrs)

    return Completion(
        operator=eff_operator, partition=part, graph=graph, exotics=exotics, terms=terms
    )


def operator_completions(operator: EffectiveOperator) -> List[Completion]:
    # print(f"Finding completions of {operator.name}...")
    parts = partitions(operator)
    completions = [partition_completion(p) for p in parts]
    good_completions = [c for c in completions if not isinstance(c, FailedCompletion)]
    return good_completions


def compare_terms(comp1: Completion, comp2: Completion) -> Dict[str, str]:
    """Returns a dictionary representing the field relabellings that would need to
    be applied to the terms of comp1 to make it equivalent to comp2. This
    includes the identity remapping. That is, if the terms of comp1 are the same
    as the terms in comp2, the function returns a dictionary like

       {"φ": "φ", "η": "η", ...}

    """
    # make sure field content is the same
    if set(comp1.exotic_info().values()) != set(comp2.exotic_info().values()):
        return {}

    # cannot be equivalent
    if len(comp1.terms) != len(comp2.terms):
        return {}

    # qnumbers -> label
    rev_map2 = {
        qnumbers: field.label for field, qnumbers in comp2.exotic_info().items()
    }
    remapping = {
        rev_map2[qnumbers]: field.label
        for field, qnumbers in comp1.exotic_info().items()
    }

    new_terms = set()
    for term in comp1.terms:
        for k, v in remapping.items():
            simple = term.simplify()
            s = str(safe_nocoeff(simple))
            s = multiple_replace(remapping, s)
            # remove generation indices in comparison
            s = re.sub(r"g[0-9]+_", "g_", s)
            # remove negative signs on indices
            s = re.sub(r"-", "", s)
            ss = tuple(sorted(s.split("*")))
            new_terms.add(ss)

    comp2_strs = []
    for term in comp2.terms:
        simple = term.simplify()
        s = str(safe_nocoeff(simple))
        comp2_strs.append(s)

    comp2_strs = [re.sub(r"g[0-9]+_", "g_", s) for s in comp2_strs]
    comp2_strs = [re.sub(r"-", "", s) for s in comp2_strs]
    comp2_tups = [tuple(sorted(s.split("*"))) for s in comp2_strs]

    if new_terms == set(comp2_tups):
        return remapping

    # otherwise, no equivalence, return empty dict
    return {}


def are_equivalent_completions(comp1: Completion, comp2: Completion) -> bool:
    """Checks to see if the Lagrangian terms describing two completions are
    equivalent.

    Two completions are equivalent if their Lagrangian terms in canonical form
    are the same up to field relabellings.

    """
    return bool(compare_terms(comp1, comp2))


def remove_equivalent_completions(comps: List[Completion]) -> List[Completion]:
    """Compares completions by comparing Lagrangian terms. Removes duplicates and
    returns copied list.

    """
    remove_equivalent(comps, are_equivalent_completions)


# TODO change_to_collect_models
def collect_completions(
    completions: List[Completion], key=None
) -> Dict[tuple, Completion]:
    """Return dictionary mapping field content to list of completions.

    `key` is a function that takes a completion and returns a dictionary mapping
    FieldType to a tuple of numbers representing that field. This defaults to
    the `exotic_info` method.

    """
    out = {}

    if key is None:
        key = lambda x: x.exotic_info()

    func = lambda c: tuple(sorted(key(c).values()))
    for k, g in groupby(completions, key=func):
        g_list = list(g)
        # TODO calling this here will be a bottleneck
        remove_equivalent_completions(g_list)
        k = tuple(sorted(set(k)))
        out[k] = g_list

    return out


def prime_registry(sieve: Dict[tuple, List[Completion]]) -> Dict[tuple, int]:
    """Ascribe a unique prime number to each exotic appearing in `sieve`.

    `sieve` is a dictionary mapping a tuple of field information to a list of
    completions.

    """
    reg = {}
    counter = 1
    for k, v in sieve.items():
        for field in k:
            if field not in reg:
                reg[field] = prime(counter)
                counter += 1
    return reg


def model_registry(completions, registry) -> Dict[tuple, int]:
    """Assigns an unique integer to every model by multiplying primes of fields."""
    reg = {}
    for k in completions:
        prod = 1
        for field in k:
            prod *= registry[field]

        reg[k] = prod

    return reg


# TODO Change to filter_models
def filter_completions(
    completions: Dict[tuple, List[Completion]], sieve: Dict[tuple, List[Completion]]
) -> Dict[tuple, List[Completion]]:
    # establish prime registry
    registry = prime_registry({**sieve, **completions})

    # construct dictionaries mapping tuples of field info to integers (products
    # of primes)
    completions_model_registry = model_registry(completions, registry)
    sieve_model_registry = model_registry(sieve, registry)

    unique = {}
    for k, v in completions_model_registry.items():
        factors_ = factors(v)
        for ref_val in sieve_model_registry.values():
            if ref_val in factors_:
                break
        else:  # no break => unique model
            unique[k] = completions[k]

    return unique