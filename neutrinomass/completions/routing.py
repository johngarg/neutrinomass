"""Independent physical endpoint validation for scalar/fermion trees.

Legacy graphs are accepted only when occurrence inference is unambiguous.
Construction records explicit ports before species bookkeeping loses identity.
"""

from collections import defaultdict
from dataclasses import dataclass

import networkx as nx

from neutrinomass.completions.amplitudes import (
    external_leg_provenance, _internal_occurrence_edges, _interaction_vertices,
)
from neutrinomass.completions.completions import canonical_propagator_cut, partition_leaves
from neutrinomass.completions.core import VectorLikeDiracFermion
from neutrinomass.completions.equivalence import field_label_parts


class InvalidMomentumGraph(ValueError):
    """The final UV graph does not support its claimed momentum routing."""


def validate_routed_completion(completion):
    """Consumer boundary: no cached certification of mutable graph metadata."""
    if completion.momentum_contributions:
        return validate_routes(completion)
    if completion.derivative_edges:
        raise InvalidMomentumGraph("Legacy edge-only routing requires explicit migration")


def state_key(field, ignored_charges=()):
    """Actual field state, without dummy indices or derivative helper labels."""
    field = field.strip_derivs() if field.derivs else field
    # strip_derivs returns an unindexed Field, whose conjugation lives in
    # is_conj rather than its label. IndexedField labels already carry it.
    label = field.label
    if field.is_conj and not label.endswith("†"):
        label += "†"
    charges = tuple(sorted(
        (name, value) for name, value in field.charges.items()
        if name not in ignored_charges
    ))
    return (label, field.dynkin, charges, field.is_conj)


def species_states(exotics, ignored_charges=()):
    """Resolve states by explicit conjugation/partner operations, not suffixes."""
    states = {}
    for exotic in exotics:
        orbit = [exotic, exotic.conj]
        if isinstance(exotic, VectorLikeDiracFermion):
            orbit.extend([exotic.dirac_partner(), exotic.dirac_partner().conj])
        family = frozenset(state_key(f, ignored_charges) for f in orbit)
        for field in orbit:
            key = state_key(field, ignored_charges)
            value = (family, type(exotic).__name__, getattr(field, "is_unbarred", None))
            if key in states and states[key] != value:
                raise InvalidMomentumGraph("Ambiguous species/partner state")
            states[key] = value
    return states


def _conjugate_gauge_states(left, right, ignored_charges):
    return left.dynkin[2:] == right.conj.dynkin[2:] and all(
        left.charges.get(k, 0) + right.charges.get(k, 0) == 0
        for k in (left.charges.keys() | right.charges.keys()) - set(ignored_charges)
    )


def classify_fermion_pair(left, right, states, ignored_charges=()):
    """Check the Majorana or vectorlike mass/kinetic pairing independently."""
    if not left.is_fermion or not right.is_fermion:
        raise InvalidMomentumGraph("A fermion propagator has a nonfermion endpoint")
    try:
        a = states[state_key(left, ignored_charges)]
        b = states[state_key(right, ignored_charges)]
    except KeyError as error:
        raise InvalidMomentumGraph("Endpoint is not an explicit registered Weyl state") from error
    if a[:2] != b[:2]:
        raise InvalidMomentumGraph("Propagator connects different species")
    if not _conjugate_gauge_states(left, right, ignored_charges):
        raise InvalidMomentumGraph("Propagator gauge representations/charges are not conjugate")
    momentum = left.lorentz_irrep != right.lorentz_irrep
    if left.dynkin[:2] not in {"10", "01"} or right.dynkin[:2] not in {"10", "01"}:
        raise InvalidMomentumGraph("Only Weyl fermion propagators are supported")
    if a[1] == "VectorLikeDiracFermion":
        if (a[2] == b[2]) != momentum:
            raise InvalidMomentumGraph("Wrong vectorlike Weyl partner for mass/kinetic propagation")
    elif a[1] != "MajoranaFermion":
        raise InvalidMomentumGraph("Unrecognised heavy fermion reality")
    return "momentum" if momentum else "mass"


@dataclass(frozen=True)
class EndpointAudit:
    occurrences: dict
    vertices: dict
    momentum_edges: tuple
    stored_momentum_edges: tuple
    edge_kinds: dict

    @property
    def agrees(self):
        return self.momentum_edges == self.stored_momentum_edges


def inspect_endpoints(completion):
    if not nx.is_tree(completion.graph):
        raise InvalidMomentumGraph("Matching requires a connected tree")
    external = partition_leaves(completion.partition)
    # 3B is a bookkeeping charge. On a B-violating tree the two sides of
    # a cut can infer different 3B assignments for the same gauge multiplet.
    total_3b = sum(leaf.field.charges.get("3b", 0) for leaf in external)
    ignored_charges = ("3b",) if total_3b else ()
    states = species_states(completion.exotics, ignored_charges)
    labels = {label for key in states for label in [key[0]]}
    bases = {field_label_parts(label)[0] for label in labels}
    provenance = external_leg_provenance(completion)
    vertices = _interaction_vertices(completion, provenance, bases)
    # Legacy label-to-port inference is valid only with one incident edge of
    # each species. Repeated-species graphs must carry explicit heavy_ports.
    for vertex in vertices.values():
        if "heavy_ports" in completion.graph.nodes[vertex]:
            continue
        incident = [field_label_parts(d["particle"])[0]
                    for _, _, d in completion.graph.edges(vertex, data=True)
                    if field_label_parts(d["particle"])[0] in bases]
        if len(incident) != len(set(incident)):
            raise InvalidMomentumGraph("Ambiguous repeated-species ports: explicit occurrence provenance required")
    occurrences = _internal_occurrence_edges(completion, provenance, bases)
    expected = {(ti, fi) for ti, term in enumerate(completion.terms)
                for fi, field in enumerate(term.tensors)
                if hasattr(field, "label") and field_label_parts(field.label)[0] in bases}
    if set(occurrences) != expected:
        raise InvalidMomentumGraph("Heavy occurrence provenance does not cover the terms exactly")
    endpoints = defaultdict(list)
    for occurrence, edge in occurrences.items():
        endpoints[tuple(sorted(edge))].append(occurrence)
    leaves = {leaf.node for leaf in external}
    internal = {tuple(sorted((a, b))) for a, b in completion.graph.edges
                if a not in leaves and b not in leaves}
    if set(endpoints) != internal:
        raise InvalidMomentumGraph("Internal edges and occurrence provenance differ")
    kinds = {}
    for edge, ports in endpoints.items():
        if len(ports) != 2 or {vertices[ti] for ti, _ in ports} != set(edge):
            raise InvalidMomentumGraph("Heavy edge does not have exactly one port at each endpoint")
        left, right = [completion.terms[ti].tensors[fi] for ti, fi in ports]
        edge_species = field_label_parts(completion.graph.edges[edge]["particle"])[0]
        if any(field_label_parts(f.label)[0] != edge_species for f in (left, right)):
            raise InvalidMomentumGraph("Graph species label disagrees with its endpoints")
        if left.derivs or right.derivs:
            raise InvalidMomentumGraph("Differentiated heavy endpoints require the signed-vertex route validator")
        if left.is_fermion or right.is_fermion:
            kinds[edge] = classify_fermion_pair(left, right, states, ignored_charges)
        else:
            try:
                a = states[state_key(left, ignored_charges)]
                b = states[state_key(right, ignored_charges)]
            except KeyError as error:
                raise InvalidMomentumGraph("Unregistered bosonic endpoint") from error
            if a[:2] != b[:2] or left.lorentz_irrep != right.lorentz_irrep:
                raise InvalidMomentumGraph("Bosonic endpoint species/spins disagree")
            if not _conjugate_gauge_states(left, right, ignored_charges):
                raise InvalidMomentumGraph("Bosonic endpoint gauge states are not conjugate")
            if a[1].startswith("Complex") and left.is_conj == right.is_conj:
                raise InvalidMomentumGraph("Complex propagator has equally oriented endpoints")
            kinds[edge] = "scalar"
    momentum = tuple(sorted(e for e, kind in kinds.items() if kind == "momentum"))
    stored = tuple(sorted(tuple(sorted(p.edge)) for p in completion.momentum_contributions
                          if p.numerator_kind == "momentum"))
    return EndpointAudit(occurrences, vertices, momentum, stored, kinds)


def validate_routes(completion, *, derivative_order=None):
    """Validate before application projection, equivalence or persistence."""
    audit = inspect_endpoints(completion)
    if not audit.agrees:
        raise InvalidMomentumGraph(f"Stored momentum edges {audit.stored_momentum_edges} disagree with final endpoints {audit.momentum_edges}")
    contributions = completion.momentum_contributions
    routes = tuple(p.derivative_route() for p in contributions if p.numerator_kind == "momentum")
    if routes != completion.derivative_routes or tuple(p.edge for p in contributions) != completion.derivative_edges:
        raise InvalidMomentumGraph("Mutable route/contribution/edge views disagree")
    if len({tuple(sorted(p.edge)) for p in contributions}) != len(contributions):
        raise InvalidMomentumGraph("Duplicate propagator contributions on an edge")
    for p in contributions:
        edge = tuple(sorted(p.edge))
        if p.numerator_kind != audit.edge_kinds.get(edge):
            raise InvalidMomentumGraph("Stored numerator kind disagrees with final endpoints")
        if p.particle != completion.graph.edges[edge]["particle"]:
            raise InvalidMomentumGraph("Stored numerator species disagrees with the verified edge")
        if p.cut_side and tuple(p.cut_side) != canonical_propagator_cut(completion.partition, completion.graph, edge):
            raise InvalidMomentumGraph("Momentum cut is not the canonical cut of the verified edge")
    vertex_degree = sum(f.derivs for term in completion.terms for f in term.indexed_fields)
    if derivative_order is not None and vertex_degree + sum(p.derivative_degree for p in contributions) != derivative_order:
        raise InvalidMomentumGraph("Vertex + numerator + denominator degree exceeds or misses the requested total budget")
    return audit


def migrate_legacy_one_derivative(completion):
    """Repair a complete, unambiguous raw graph; invalidate old projection.

    This cannot certify coverage of a stream pruned with the old polynomial.
    The caller must retain the old bytes, hash and producer separately.
    """
    from copy import copy, deepcopy
    audit = inspect_endpoints(completion)
    contributions = completion.momentum_contributions
    if (len(audit.momentum_edges) != 1 or len(contributions) != 1
            or contributions[0].numerator_kind != "momentum"
            or contributions[0].denominator_order != 0
            or any(f.derivs for t in completion.terms for f in t.indexed_fields)):
        raise InvalidMomentumGraph("Migration requires one kinetic numerator and no derivative vertices")
    migrated = copy(completion)
    migrated.graph = deepcopy(completion.graph)
    edge = audit.momentum_edges[0]
    particle = completion.graph.edges[edge]["particle"]
    port = next(o for o, e in audit.occurrences.items() if tuple(sorted(e)) == edge)
    field = completion.terms[port[0]].tensors[port[1]]
    contribution = contributions[0]._replace(
        edge=edge, particle=particle, numerator_lorentz=field.dynkin[:2],
        cut_side=canonical_propagator_cut(completion.partition, completion.graph, edge),
    )
    migrated.momentum_contributions = (contribution,)
    migrated.derivative_routes = (contribution.derivative_route(),)
    migrated.derivative_edges = (edge,)
    migrated.lorentz_projection = None
    for ti, vertex in audit.vertices.items():
        migrated.graph.nodes[vertex]["interaction"] = ti
        migrated.graph.nodes[vertex]["heavy_ports"] = [
            {"tensor": fi, "edge": list(e)}
            for (term, fi), e in sorted(audit.occurrences.items()) if term == ti
        ]
    validate_routes(migrated, derivative_order=1)
    return migrated

