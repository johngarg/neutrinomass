"""Construct local propagator expansions from the finished physical UV tree."""

from itertools import product
import networkx as nx
from sympy import Matrix, Rational

from neutrinomass.completions.core import PropagatorContribution
from neutrinomass.completions.completions import (
    _basis_singlets, _completion_with_momentum_contributions,
    _edge_particle_field, _embed_multi_derivative_contraction,
    _stripped_occurrence_key, canonical_propagator_cut,
    expand_propagator_denominators, partition_leaves,
)
from neutrinomass.completions.routing import inspect_endpoints, validate_routes
from neutrinomass.tensormethod.core import Operator, eps
from neutrinomass.tensormethod.lorentz import lorentz_contraction, primitive_coordinates


def physical_numerators(completion):
    audit = inspect_endpoints(completion)
    return tuple(PropagatorContribution(
        edge=edge, particle=_edge_particle_field(completion,edge)[0],
        numerator_kind="momentum", denominator_order=0,
        numerator_lorentz=_edge_particle_field(completion,edge)[1].dynkin[:2],
        differentiated_field="", differentiated_lorentz="",
        cut_side=canonical_propagator_cut(completion.partition,completion.graph,edge))
        for edge in audit.momentum_edges)


def sewn_fermion_operator(completion, momentum_fields=()):
    """Sew mass chains and expose a labelled bispinor for each kinetic edge.

    Renormalizable scalar/fermion vertices contain zero or two Weyl ports.
    Removing kinetic edges leaves paths pairing external and spurion ports;
    a path with incompatible chiralities is a physical construction error.
    Overall epsilon orientation is fixed deterministically; projection
    conventions identify a common nonzero normalization.
    """
    audit = inspect_endpoints(completion)
    if len(momentum_fields) != len(audit.momentum_edges):
        raise ValueError("A momentum spurion is required for every kinetic edge")
    external = {leaf.node: leaf.field.strip_derivs_with_indices()
                for leaf in partition_leaves(completion.partition)}
    chains = nx.Graph()
    ends = {}
    for node, field in external.items():
        if field.is_fermion:
            chains.add_edges_from(completion.graph.edges(node))
            kind = "Undotted" if field.dynkin[:2] == "10" else "Dotted"
            index, = field.indices_by_type[kind]
            ends[node] = (kind,index)
    chains.add_edges_from(edge for edge,kind in audit.edge_kinds.items() if kind == "mass")
    momenta = []
    occupied={str(index) for field in external.values() for index in field.indices}
    for ordinal,(edge, field) in enumerate(zip(audit.momentum_edges,momentum_fields)):
        # Loaded records do not advance the process-global fresh-index
        # counters. Use explicitly reserved, checked names rather than risk
        # colliding a spurion with an archived external spinor index.
        suffix=f"physical_momentum_{ordinal}"
        while {"u"+suffix,"d"+suffix}&occupied:
            suffix += "x"
        momentum = field("u"+suffix+" d"+suffix)
        occupied.update(str(index) for index in momentum.indices)
        momenta.append(momentum)
        for vertex in edge:
            ports = [o for o,e in audit.occurrences.items()
                     if tuple(sorted(e)) == edge and audit.vertices[o[0]] == vertex]
            if len(ports) != 1:
                raise ValueError("Kinetic endpoint has ambiguous tensor provenance")
            ti,fi = ports[0]
            endpoint = completion.terms[ti].tensors[fi]
            kind = "Undotted" if endpoint.dynkin[:2] == "10" else "Dotted"
            index, = momentum.indices_by_type[kind]
            key = ("momentum",edge,vertex)
            ends[key] = (kind,index)
            chains.add_edge(vertex,key)
    contractions = []
    for component in nx.connected_components(chains):
        terminals = sorted(component & ends.keys(),key=repr)
        if len(terminals) != 2 or any(chains.degree(v) > 2 for v in component):
            raise ValueError("UV Weyl chains must pair exactly two exposed spinor ports")
        (left,a),(right,b) = [ends[node] for node in terminals]
        if left != right:
            raise ValueError("A mass/Yukawa chain joins incompatible chiralities")
        contractions.append(eps(f"{-a} {-b}"))
    return Operator(*external.values(),*momenta,*contractions)


def _box_projection(basis, contraction, placement):
    """Project the scalar D^2 component by tracing both momentum spinors.

    The embedded box columns span the antisymmetric momentum pair in each
    SL(2) factor. Their Gram projection removes the symmetric, traceless
    derivative irreps without imposing an equation of motion.
    """
    assignments = {placement.occurrence_keys[0]: tuple(f.label for f in basis.momentum_fields)}
    columns = [_embed_multi_derivative_contraction(
        lorentz_contraction(singlet,normalise=True), placement.operator.operator,
        basis.ambient_operator, assignments).evaluation_vector()
        for singlet in _basis_singlets(placement.operator.operator,placement.basis)]
    matrix = Matrix.hstack(*(Matrix(column) for column in columns))
    gram = matrix.T*matrix
    return tuple(gram.inv()*matrix.T*Matrix(contraction.evaluation_vector()))


def project_second_derivative_tree(basis, completion):
    """Project a physical p^2 denominator or two sewn kinetic numerators."""
    validate_routes(completion,derivative_order=2)
    momenta = [p for p in completion.momentum_contributions if p.numerator_kind == "momentum"]
    if not momenta:
        denominator, = completion.momentum_contributions
        return basis.project_denominator(sewn_fermion_operator(completion),
                                         completion.partition,denominator)
    if len(momenta) != 2 or any(p.denominator_order for p in momenta):
        raise ValueError("Unexpected order-two propagator expansion")
    audit = inspect_endpoints(completion)
    by_edge = {tuple(sorted(p.edge)):p for p in momenta}
    cuts = [by_edge[edge].cut_side for edge in audit.momentum_edges]
    contraction = lorentz_contraction(sewn_fermion_operator(completion,basis.momentum_fields))
    external = {leaf.node:basis.occurrence_keys.index(_stripped_occurrence_key(leaf.field))
                for leaf in partition_leaves(completion.partition)}
    coordinates = [Rational(0)]*len(basis.labels)
    for nodes in product(*cuts):
        indices = tuple(external[node] for node in nodes)
        placement = basis.placement_by_indices.get(tuple(sorted(indices)))
        if placement is None:
            continue
        local = (_box_projection(basis,contraction,placement) if indices[0] == indices[1]
                 else basis._project_cross(contraction,indices))
        for offset,value in enumerate(local):
            coordinates[placement.coordinate_start+offset] += value
    coordinates = primitive_coordinates(coordinates)
    return basis._metadata(coordinates) if any(coordinates) else None


def local_physical_expansions(completion, derivative_order, *, second_derivative_basis=None):
    numerators = physical_numerators(completion)
    remaining = derivative_order-len(numerators)
    if remaining < 0 or remaining % 2:
        return
    base = _completion_with_momentum_contributions(completion,numerators)
    for candidate in expand_propagator_denominators(base,remaining//2):
        validate_routes(candidate,derivative_order=derivative_order)
        if second_derivative_basis is not None:
            projection = project_second_derivative_tree(second_derivative_basis,candidate)
            if projection is None:
                continue
            candidate.lorentz_projection = projection
        yield candidate
