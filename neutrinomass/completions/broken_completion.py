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
    return


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


def lorentz_pairs_maybe_remove(
    left: IndexedField, right: IndexedField, free_indices: List[Index]
):
    """Returns a list of epsilons of paired dotted and undotted indices, and removes
    these from free_indices by side-effect if free_indices is not empty.

    """

    lu, ld, _, _, _ = left.indices_by_type.values()
    ru, rd, _, _, _ = right.indices_by_type.values()

    out = []
    index_pairs = list(zip(lu, ru)) + list(zip(ld, rd))
    for i, j in index_pairs:
        epsilon = eps(f"-{i} -{j}")
        out.append(epsilon)
        if free_indices:
            free_indices.remove(i)
            free_indices.remove(j)

    return out


# def scalar_contraction_epsilons(
#     undotted, dotted, two_undotted, two_dotted, free_indices
# ):
#     lorentz_contraction_epsilons = []
#     to_remove = []
#     if two_undotted:
#         eps_indices = " ".join(str(-i) for i in undotted)
#         lorentz_eps = eps(eps_indices)
#         lorentz_contraction_epsilons.append(lorentz_eps)
#         to_remove += undotted

#     if two_dotted:
#         eps_indices = " ".join(str(-i) for i in dotted)
#         lorentz_eps = eps(eps_indices)
#         lorentz_contraction_epsilons.append(lorentz_eps)
#         to_remove += dotted

#     if not (two_undotted or two_dotted):
#         return []

#     if free_indices:
#         for i in to_remove:
#             free_indices.remove(i)

#     return lorentz_contraction_epsilons


def get_coupling_field(
    fields: Tuple[IndexedField],
    symbols: Dict[str, List[str]],
    field_dict: Dict[tuple, str],
    n_derivs: int,
    free_indices: List[Index],
) -> Tuple[FieldType, FieldType]:
    """Return the new exotic field that couples to `fields` along with its partner
    that features in the Lagrangian term.

    free_indices passed in separately because it has been mutated to remove the
    contracted Lorentz indices.

    """

    # construct exotic field
    # deal with charges
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
    if lorentz_dynkin in ("10", "01"):
        symbols_to_use = symbols["fermion"]
    else:
        symbols_to_use = symbols["boson"]

    field_labels = sorted([f.field for f in fields], key=lambda x: x.label_with_dagger)
    field_labels = tuple(field_labels)
    if field_labels in field_dict.keys():
        symbol = field_dict[field_labels]
    else:
        symbol = symbols_to_use.pop(0)
        field_dict[field_labels] = symbol

    # for Dirac and Majorana fermions, always keep plain symbol left handed
    to_conj = False
    if exotic_indices and exotic_indices[0] == "d":
        exotic_indices = Index.conj_index_string(exotic_indices)
        to_conj = True

    exotic_field = IndexedField(
        label=symbol, indices=exotic_indices, charges=exotic_charges
    )

    # construct MajoranaFermion, VectorLikeDiracFermion, ...
    # `partner` is in the Lagrangian, `exotic_field` transforms like contracted pair
    exotic_field = cons_completion_field(exotic_field)
    exotic_field = exotic_field.conj if to_conj else exotic_field

    if isinstance(exotic_field, VectorLikeDiracFermion):
        partner = exotic_field.dirac_partner() if not n_derivs else exotic_field.conj
    elif isinstance(exotic_field, ComplexScalar):
        partner = exotic_field.conj
    elif isinstance(exotic_field, MajoranaFermion):
        partner = exotic_field.majorana_partner()
    else:
        raise Exception(f"Unrecognised exotic field type: {type(exotic_field)}")

    return exotic_field, partner


def is_allowed_lorentz_structure(
    dotted: List[Index], undotted: List[Index], n_derivs: int
) -> Tuple[bool, bool, bool]:
    # sort out lorentz epsilons from contraction
    # contractions between scalars
    scalar = len(undotted) == 0 and len(dotted) == 0
    # contractions between undotted fermions (or undotted and one-deriv dotted)
    two_undotted = len(undotted) == 2 and len(dotted) < 2
    # contractions between dotted fermions (or dotted and one-deriv undotted)
    two_dotted = len(dotted) == 2 and len(undotted) < 2
    # contraction between scalar and undotted fermion (or one-deriv dotted)
    one_undotted = len(undotted) == 1 and len(dotted) == 0
    # contraction between scalar and dotted fermion (or one-deriv undotted fermion)
    one_dotted = len(dotted) == 1 and len(undotted) == 0
    # fermion derivative cases
    # derivative acting on scalar into an undotted fermion (or one-deriv dotted)
    # two_undotted_one_dotted = len(undotted) == 2 and len(dotted) == 1
    # derivative acting on scalar into a dotted fermion (or one-deriv undotted)
    # two_dotted_one_undotted = len(dotted) == 2 and len(undotted) == 1

    # Case that always simplifies two derivatives on scalar down. Corresponds to
    # the case of D...D(scalar) contracted with D...D(scalar). Completion should
    # behave the same as `scalar`
    # even_undotted_even_dotted = (
    #     len(dotted) != 0
    #     and len(undotted) != 0
    #     and len(undotted) % 2 == 0
    #     and len(dotted) % 2 == 0
    # )

    # if not an allowed contraction, return none
    if (
        scalar
        or two_undotted
        or two_dotted
        or one_undotted
        or one_dotted
        # or two_undotted_one_dotted
        # or two_dotted_one_undotted
        # or even_undotted_even_dotted
    ):
        return True, two_undotted, two_dotted

    return False, two_undotted, two_dotted


def get_deriv_term(
    fields: List[IndexedField],
    epsilons: List[Tensor],
    n_derivs: int,
    extra_indices: List[Index],
) -> Operator:
    """Takes a list of fields with free Lorentz indices and constructs Lagrangian
    term from which the derivative contraction may have arisen. (This should be
    unique.)

    """

    if len(extra_indices) == 2:
        epsilons.append(eps(" ".join(str(-i) for i in extra_indices)))

    # Corresponds to cases (DF)(DS) and (DS)(DS), (DF)(DF)
    no_deriv_fields = []
    for f in fields:
        # if field doesn't have a derivative, just add it to the list
        if f.stripped is None:
            no_deriv_fields.append(f)
            continue

        head = f.strip_derivs()
        gauge_indices = f.gauge_indices

        if f.is_fermion:
            fermion_index = extra_indices.pop(0)
            if n_derivs == 1:
                fermion_index = fermion_index.conj
            indices = [fermion_index, *gauge_indices]
        else:
            indices = gauge_indices

        no_deriv_indexed_field = head(" ".join(str(i) for i in indices))
        no_deriv_fields.append(no_deriv_indexed_field)

    if n_derivs == 2:
        return no_deriv_term(no_deriv_fields, epsilons)

    # Corresponds to cases (DS)f, (Df)S
    # breakpoint()

    # construct term
    new_epsilons = [*contracted_epsilons, *lorentz_contraction_epsilons]

    term = reduce(lambda a, b: a * b, new_epsilons, prod_fields)
    term *= partner

    for free in term.free_indices:
        assert free.index_type == "Generation"

    return term, spectator_epsilons, lorentz_contraction_epsilons


def no_deriv_term(fields, epsilons):
    # construct term
    prod_fields = reduce(lambda x, y: x * y, fields)
    term = reduce(lambda a, b: a * b, epsilons, prod_fields)
    return term


def contract(
    fields: Tuple[IndexedField],
    symbols: Dict[str, List[str]],
    epsilons: List[Tensor],
    field_dict: Dict[tuple, str],
) -> Tuple[FieldType, Operator, List[Tensor], List[Tensor]]:
    """Takes two or three indexed fields and the epsilons and returns (amongst other
    things) a new indexed field transforming in the same way as $x otimes y$.

    Returns a possibly empty tuple of: the new field, the Lagrangian term, the
    epsilons passed around (and possibly removed), new epsilons from contraction
    of Lorentz indices.

    Example:
        >>> field, term, epsilons, new_epsilons = contract((H('i0'), H('i1')), {"fermion": [], "boson": ["S"]}, [], {})
        >>> field
        S(i0, i1)
        >>> field.y
        1

    """
    if len(fields) != 2 and len(fields) != 3:
        raise Exception("Too many fields passed to contract.")

    # product of fields
    prod_fields = reduce(lambda a, b: a * b, fields)
    free_indices = prod_fields.get_free_indices()

    contracted_epsilons, spectator_epsilons = [], []
    # contracted epsilons are those all of whose indices are contracted on the
    # set of fields passed in. Spectator epsilons are all of the others
    for epsilon in epsilons:
        if is_contracted_epsilon(epsilon, free_indices):
            contracted_epsilons.append(epsilon)
        else:
            spectator_epsilons.append(epsilon)

    n_derivs = 0
    for f in fields:
        if str(f).startswith("D"):
            n_derivs += 1

    # scalar_contraction_epsilons removes pairs of dotted or undotted indices
    # from free_indices so that the exotic field can be found
    #
    # When only dealing with two derivatives, there are at most only two fields
    # that carry Lorentz indices, call left and right
    nontrivial_lorentz = [f for f in fields if f.lorentz_irrep != (0, 0)]
    assert len(nontrivial_lorentz) < 3
    if len(nontrivial_lorentz) < 1:
        lorentz_contraction_epsilons = []
    elif len(nontrivial_lorentz) == 1:
        # need to contract on same term, e.g. D(F) S
        for k, g in groupby(free_indices, key=lambda i: i.index_type):
            grp = list(g)
            if len(grp) == 2:
                lorentz_contraction_epsilons = [eps(f"-{grp[0]} -{grp[1]}")]
                free_indices.remove(grp[0])
                free_indices.remove(grp[1])
    elif len(nontrivial_lorentz) == 2:
        left, right = nontrivial_lorentz
        lorentz_contraction_epsilons = lorentz_pairs_maybe_remove(
            left, right, free_indices
        )
    else:
        raise Exception(
            "Don't currently support more than 2 derivatives in an operator."
        )

    undotted, dotted, _, _, _ = Index.indices_by_type(free_indices).values()
    is_allowed = is_allowed_lorentz_structure(dotted, undotted, n_derivs)
    if not is_allowed:
        return None

    exotic_field, partner = get_coupling_field(
        fields, symbols, field_dict, n_derivs, free_indices
    )

    # Need additional su2 epsilons to fix su2 indices (since not working with
    # lowered indices at all). Won't need to do this if removing a derivative in
    # the process

    # if n_derivs == 2 and exotic_field.is_fermion:
    # if there is an odd number of derivatives acting on a fermion
    partner, fix_su2_epsilons = exotic_field.lower_su2()
    # elif n_derivs == 1:
    #     partner, fix_su2_epsilons = partner.lower_su2()

    # In D(f)D(f) case, extra indices are the indices of the fermions.
    # Specialising to the case of only two derivatives in the operator, these
    # indices come from the second element of the lorentz_contraction_epsilons
    # if n_derivs:
    #     if exotic_field.is_fermion:
    #         extra_indices = [
    #             i for i in free_indices if i.index_type in ("Dotted", "Undotted")
    #         ]
    #     else:
    #         extra_indices = [-i for i in lorentz_contraction_epsilons[1].indices]

    #     term = get_deriv_term(
    #         fields + (partner,),
    #         [*lorentz_contraction_epsilons, *contracted_epsilons, *fix_su2_epsilons],
    #         n_derivs,
    #         extra_indices,
    #     )
    # else:
    new_epsilons = [
        *contracted_epsilons,
        *fix_su2_epsilons,
        *lorentz_contraction_epsilons,
    ]
    # breakpoint()
    term = no_deriv_term(fields + (partner,), new_epsilons)

    for free in term.free_indices:
        assert free.index_type == "Generation"

    return term, spectator_epsilons, lorentz_contraction_epsilons


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
    """

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

    # contractions between undotted fermions (or undotted and one-deriv dotted)
    two_undotted = len(undotted) == 2 and len(dotted) == 0
    # contractions between dotted fermions (or dotted and one-deriv undotted)
    two_dotted = len(dotted) == 2 and len(undotted) == 0
    # fermion derivative cases
    # derivative acting on scalar into an undotted fermion (or one-deriv dotted)
    two_undotted_one_dotted = len(undotted) == 2 and len(dotted) == 1
    # derivative acting on scalar into a dotted fermion (or one-deriv undotted)
    two_dotted_one_undotted = len(dotted) == 2 and len(undotted) == 1
    # Case that always simplifies two derivatives on scalar down. Corresponds to
    # the case of D...D(scalar) contracted with D...D(scalar). Completion should
    # behave the same as `scalar`
    even_undotted_even_dotted = (
        len(dotted) != 0
        and len(undotted) != 0
        and len(undotted) % 2 == 0
        and len(dotted) % 2 == 0
    )

    lorentz_contraction_epsilons = scalar_contraction_epsilons(
        undotted=undotted,
        dotted=dotted,
        two_undotted=two_undotted,
        two_dotted=two_dotted,
        free_indices=[],
    )
    for e in lorentz_contraction_epsilons:
        prod *= e

    free_indices = prod.free_indices
    contracted_epsilons, spectator_epsilons = [], []
    for e in epsilons:
        if is_contracted_epsilon(e, free_indices):
            contracted_epsilons.append(e)
        else:
            spectator_epsilons.append(e)

    prod = reduce(lambda a, b: a * b, contracted_epsilons, prod)

    for free in prod.free_indices:
        assert free.index_type == "Generation"

    terms.append(prod)

    return terms, edge_dict, field_dict, lorentz_epsilons


def partition_completion(partition) -> Union[Completion, FailedCompletion]:
    part = partition["partition"]
    epsilons = partition["epsilons"]
    graph = partition["graph"]
    op = partition["operator"]

    args = cons_completion(partition=part, epsilons=epsilons, graph=graph)
    if args is not None:
        terms, edge_dict, field_dict, lorentz_epsilons = args
    else:
        return FailedCompletion("Bad lorentz irrep")

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

    # pmatch_fields = [f.pmatch_data for f in fields]

    # # sort fields so that boson comes before fermion
    # fields = sorted(fields, key=lambda f: f.comm)

    # # Note: two boson or two fermions are not allowed lorentz structures in
    # # a contraction with one derivative

    # # match on (DS)f case
    # ds_f = pmatch(
    #     [
    #         (
    #             "?ss",
    #             (
    #                 ("Undotted", "?us"),
    #                 ("Dotted", "?ds"),
    #                 ("Colour", "?cs"),
    #                 ("Isospin", "?is"),
    #                 ("Generation", "?gs"),
    #             ),
    #             (("3b", "?3bs"), ("y", "?ys")),
    #             "bose",
    #             1,
    #             "?conjs",
    #         ),
    #         (
    #             "?sf",
    #             (
    #                 ("Undotted", "?uf"),
    #                 ("Dotted", "?df"),
    #                 ("Colour", "?cf"),
    #                 ("Isospin", "?if"),
    #                 ("Generation", "?gf"),
    #             ),
    #             (("3b", "?3bf"), ("y", "?yf")),
    #             "fermi",
    #             0,
    #             "?conjf",
    #         ),
    #     ],
    #     pmatch_fields,
    # )

    # # match on f(DS) case
    # s_df = pmatch(
    #     [
    #         (
    #             "?ss",
    #             (
    #                 ("Undotted", "?us"),
    #                 ("Dotted", "?ds"),
    #                 ("Colour", "?cs"),
    #                 ("Isospin", "?is"),
    #                 ("Generation", "?gs"),
    #             ),
    #             (("3b", "?3bs"), ("y", "?ys")),
    #             "bose",
    #             0,
    #             "?conjs",
    #         ),
    #         (
    #             "?sf",
    #             (
    #                 ("Undotted", "?uf"),
    #                 ("Dotted", "?df"),
    #                 ("Colour", "?cf"),
    #                 ("Isospin", "?if"),
    #                 ("Generation", "?gf"),
    #             ),
    #             (("3b", "?3bf"), ("y", "?yf")),
    #             "fermi",
    #             1,
    #             "?conjf",
    #         ),
    #     ],
    #     pmatch_fields,
    # )

    # # sanity
    # assert not (ds_f and s_df)

    # if ds_f:
    #     to_skip = ["Dotted", "Undotted"]
    #     if not ds_f["?conjf"]:
    #         breakpoint()
    #         fermion_index = str(ds_f["?uf"][0])
    #         partner_index = str(partner.indices_by_type["Undotted"][0])
    #     else:
    #         fermion_index = str(ds_f["?df"][0])
    #         partner_index = str(partner.indices_by_type["Dotted"][0])

    #     no_deriv_scalar = fields[0].strip_derivs()
    #     scalar_gauge_indices = " ".join(str(i) for i in ds_f["?cs"] + ds_f["?is"])
    #     fields[0] = no_deriv_scalar(scalar_gauge_indices)
    #     lorentz_contraction_epsilons.append(
    #         eps(f"-{fermion_index} -{partner_index}")
    #     )

    # elif s_df:
    #     if not s_df["?conjf"]:
    #         partner_index = str(partner.indices_by_type["Undotted"][0])
    #     else:
    #         partner_index = str(partner.indices_by_type["Dotted"][0])

    #     no_deriv_fermion = fields[1].strip_derivs()
    #     fermion_gauge_indices = " ".join(str(i) for i in s_df["?cf"] + s_df["?if"])
    #     fields[1] = no_deriv_fermion(partner_index + " " + fermion_gauge_indices)
    #     breakpoint()

    # prod_fields = reduce(lambda a, b: a * b, fields)
