#!/usr/bin/env python3

"""Exact gauge-tensor audit of tree-level completion amplitudes.

The completion builder validates each UV interaction separately.  This module
performs the complementary check after the heavy fields have been contracted:
it asks whether Bose/Fermi symmetrisation of the complete induced tensor makes
the amplitude vanish.  Derivative placements are retained as distinct formal
field labels, so identical undifferentiated Higgs factors are symmetrised while
``H`` and ``DH`` are not accidentally identified.
"""

from collections import Counter, defaultdict, deque
from math import factorial
from typing import Dict, NamedTuple, Tuple

from neutrinomass.completions.completions import (
    canonical_propagator_cut,
    exotic_species,
    partition_leaves,
    weak_compositions,
)
from neutrinomass.completions.equivalence import field_label_parts
from neutrinomass.tensormethod.core import Index, IndexedField, Operator, delta


GAUGE_INDEX_TYPES = {"Colour", "Isospin", "Generation"}
LORENTZ_INDEX_TYPES = {"Undotted", "Dotted"}


class UnsupportedAmplitudeAudit(ValueError):
    """The stored completion lacks enough provenance for an exact audit."""


class ExternalLegProvenance(NamedTuple):
    """Stable location of one external field occurrence in the UV terms."""

    node: int
    term_index: int
    tensor_index: int
    field_label: str


class AmplitudeAudit(NamedTuple):
    """Result of the full induced-tensor symmetrisation check."""

    status: str
    reason: str
    derivative_sectors: int
    witness: str

    @property
    def is_zero(self):
        return self.status == "zero"


def _undifferentiated_label(field):
    if not field.derivs:
        return field.label
    label = field.stripped["label"]
    if field.is_conj and not label.endswith("†"):
        label += "†"
    return label


def _occurrence_key(field, *, strip_derivatives=False):
    label = (
        _undifferentiated_label(field)
        if strip_derivatives
        else field.label
    )
    return (
        label,
        tuple(
            str(index)
            for index in field.indices
            if index.index_type in GAUGE_INDEX_TYPES
        ),
    )


def external_leg_provenance(completion):
    """Map every external UV field factor to its graph leaf.

    Gauge and generation indices are inherited from the effective-operator
    partition and are stable across JSON serialisation.  A queue handles the
    rare case of several index-free fields with the same label without relying
    on Python object identity.
    """

    leaves = defaultdict(deque)
    sorted_leaves = sorted(
        partition_leaves(completion.partition), key=lambda item: item.node
    )
    for leaf in sorted_leaves:
        leaves[_occurrence_key(leaf.field, strip_derivatives=True)].append(leaf)

    exotic_labels = set(exotic_species(completion))
    provenance = {}
    for term_index, term in enumerate(completion.terms):
        for tensor_index, tensor in enumerate(term.tensors):
            if not isinstance(tensor, IndexedField):
                continue
            if field_label_parts(tensor.label)[0] in exotic_labels:
                continue
            key = _occurrence_key(tensor)
            if not leaves[key]:
                raise UnsupportedAmplitudeAudit(
                    "Cannot match an external UV field to a partition leaf: "
                    f"term {term_index}, factor {tensor_index}, {tensor}"
                )
            leaf = leaves[key].popleft()
            provenance[(term_index, tensor_index)] = ExternalLegProvenance(
                node=leaf.node,
                term_index=term_index,
                tensor_index=tensor_index,
                field_label=leaf.field.label,
            )

    unmatched = [leaf for queue in leaves.values() for leaf in queue]
    if unmatched:
        raise UnsupportedAmplitudeAudit(
            f"Unmatched external partition leaves: {unmatched}"
        )
    return provenance


def _freshened_tensors(term):
    """Freshen term-local dummy indices before multiplying UV vertices."""

    replacements = []
    seen = set()
    for tensor in term.tensors:
        for index in tensor.indices:
            if index.index_type == "Generation":
                continue
            positive = index if index.is_up else -index
            key = positive.index_type, positive.label
            if key in seen:
                continue
            seen.add(key)
            fresh = Index.fresh(Index.get_index_labels()[positive.index_type])
            replacements.extend(((positive, fresh), (-positive, -fresh)))
    return tuple(tensor.fun_eval(*replacements) for tensor in term.tensors)


def _decorated_external_field(field, derivative_degree):
    indices = [
        index for index in field.indices if index.index_type in GAUGE_INDEX_TYPES
    ]
    # SymPy 1.2 cannot construct a tensor head whose only slot is a generation
    # index (its empty gauge-symmetry direct product raises IndexError).  Such a
    # factor is a spectator to the gauge-index identity being tested here.
    if not any(index.index_type in {"Colour", "Isospin"} for index in indices):
        return None
    prefix = "" if not derivative_degree else (
        "D" if derivative_degree == 1 else f"D{derivative_degree}"
    )
    return IndexedField(
        label=prefix + field.label,
        indices=" ".join(map(str, indices)),
        charges=field.charges,
        is_conj=field.is_conj,
        symmetry=None,
        comm=field.comm,
        latex=field.latex,
        nf=field.nf,
        derivs=0,
    )


def _propagator_connectors(left, right):
    connectors = []
    left_by_type = Index.indices_by_type(left.indices)
    right_by_type = Index.indices_by_type(right.indices)
    for index_type in ("Colour", "Isospin"):
        left_indices = left_by_type[index_type]
        right_indices = right_by_type[index_type]
        if len(left_indices) != len(right_indices):
            raise UnsupportedAmplitudeAudit(
                f"Mismatched {index_type} propagator slots for {left} and {right}"
            )
        for left_index, right_index in zip(left_indices, right_indices):
            if index_type == "Isospin":
                connectors.append(
                    left_index.tensor_index_type.metric(left_index, right_index)
                )
            else:
                if left_index.is_up == right_index.is_up:
                    raise UnsupportedAmplitudeAudit(
                        f"Mismatched colour orientations for {left} and {right}"
                    )
                upper = left_index if left_index.is_up else right_index
                lower = right_index if left_index.is_up else left_index
                connectors.append(delta(f"{upper} {lower}"))
    return connectors


def _interaction_vertices(completion, provenance, exotic_labels):
    """Match stored UV terms to graph vertices using external-leg provenance."""

    leaf_nodes = {leaf.node for leaf in partition_leaves(completion.partition)}
    internal_vertices = set(completion.graph) - leaf_nodes
    term_vertices = {}
    external_by_term = defaultdict(list)
    for leg in provenance.values():
        external_by_term[leg.term_index].append(leg.node)

    for term_index, nodes in external_by_term.items():
        vertices = {
            next(iter(completion.graph.neighbors(node))) for node in nodes
        }
        if len(vertices) != 1:
            raise UnsupportedAmplitudeAudit(
                f"External legs of term {term_index} do not meet at one vertex"
            )
        term_vertices[term_index] = vertices.pop()

    if len(set(term_vertices.values())) != len(term_vertices):
        raise UnsupportedAmplitudeAudit(
            "Several UV terms were matched to the same interaction vertex"
        )

    remaining_terms = [
        index for index in range(len(completion.terms)) if index not in term_vertices
    ]
    remaining_vertices = sorted(internal_vertices - set(term_vertices.values()))

    def term_signature(term_index):
        return Counter(
            field_label_parts(field.label)[0]
            for field in completion.terms[term_index].indexed_fields
            if field_label_parts(field.label)[0] in exotic_labels
        )

    def vertex_signature(vertex):
        return Counter(
            field_label_parts(data["particle"])[0]
            for _, _, data in completion.graph.edges(vertex, data=True)
            if field_label_parts(data["particle"])[0] in exotic_labels
        )

    candidates = {
        term_index: [
            vertex
            for vertex in remaining_vertices
            if term_signature(term_index) == vertex_signature(vertex)
        ]
        for term_index in remaining_terms
    }

    def assign(pending, available):
        if not pending:
            return {}
        term_index = min(pending, key=lambda item: (len(candidates[item]), item))
        for vertex in candidates[term_index]:
            if vertex not in available:
                continue
            rest = assign(
                [item for item in pending if item != term_index],
                available - {vertex},
            )
            if rest is not None:
                return {term_index: vertex, **rest}
        return None

    inferred = assign(remaining_terms, set(remaining_vertices))
    if inferred is None:
        raise UnsupportedAmplitudeAudit(
            "Cannot match internal-only UV terms to interaction vertices"
        )
    term_vertices.update(inferred)
    if set(term_vertices.values()) != internal_vertices:
        raise UnsupportedAmplitudeAudit(
            "UV terms do not cover every interaction vertex"
        )
    return term_vertices


def _internal_occurrence_edges(completion, provenance, exotic_labels):
    """Map every heavy field factor to its exact propagator graph edge."""

    term_vertices = _interaction_vertices(completion, provenance, exotic_labels)
    occurrence_edges = {}
    for term_index, term in enumerate(completion.terms):
        vertex = term_vertices[term_index]
        occurrences = defaultdict(list)
        for tensor_index, tensor in enumerate(term.tensors):
            if not isinstance(tensor, IndexedField):
                continue
            label = field_label_parts(tensor.label)[0]
            if label in exotic_labels:
                occurrences[label].append(tensor_index)

        graph_edges = defaultdict(list)
        for neighbour in completion.graph.neighbors(vertex):
            data = completion.graph.edges[vertex, neighbour]
            label = field_label_parts(data["particle"])[0]
            if label in exotic_labels:
                graph_edges[label].append(tuple(sorted((vertex, neighbour))))

        if set(occurrences) != set(graph_edges):
            raise UnsupportedAmplitudeAudit(
                f"Heavy fields of term {term_index} do not match its graph vertex"
            )
        for label in sorted(occurrences):
            tensor_indices = sorted(occurrences[label])
            edges = sorted(graph_edges[label])
            if len(tensor_indices) != len(edges):
                raise UnsupportedAmplitudeAudit(
                    "Heavy-field multiplicity mismatch for "
                    f"{label} in term {term_index}"
                )
            for tensor_index, edge in zip(tensor_indices, edges):
                occurrence_edges[(term_index, tensor_index)] = edge
    return occurrence_edges


def induced_amplitude_operator(completion, derivative_degrees=None):
    """Contract all UV gauge tensors for one external derivative placement."""

    derivative_degrees = derivative_degrees or {}
    provenance = external_leg_provenance(completion)
    exotic_labels = set(exotic_species(completion))
    internal_edges = _internal_occurrence_edges(
        completion, provenance, exotic_labels
    )
    external = []
    invariants = []
    internal = defaultdict(list)

    for term_index, term in enumerate(completion.terms):
        freshened = _freshened_tensors(term)
        for tensor_index, (original, tensor) in enumerate(
            zip(term.tensors, freshened)
        ):
            if isinstance(tensor, IndexedField):
                base_label = field_label_parts(tensor.label)[0]
                if base_label in exotic_labels:
                    internal[internal_edges[(term_index, tensor_index)]].append(tensor)
                    continue
                leg = provenance[(term_index, tensor_index)]
                projected = _decorated_external_field(
                    tensor, derivative_degrees.get(leg.node, 0)
                )
                if projected is not None:
                    external.append(projected)
                continue

            index_types = {index.index_type for index in tensor.indices}
            if not index_types & LORENTZ_INDEX_TYPES:
                invariants.append(tensor)

    connectors = []
    for edge, occurrences in sorted(internal.items()):
        if len(occurrences) != 2:
            raise UnsupportedAmplitudeAudit(
                f"Propagator edge {edge} has {len(occurrences)} heavy-field slots"
            )
        connectors.extend(_propagator_connectors(*occurrences))

    return Operator(*external, *invariants, *connectors)


def _multinomial(parts):
    coefficient = factorial(sum(parts))
    for part in parts:
        coefficient //= factorial(part)
    return coefficient


def _operator_derivative_placement(completion):
    leaves = defaultdict(deque)
    sorted_leaves = sorted(
        partition_leaves(completion.partition), key=lambda item: item.node
    )
    for leaf in sorted_leaves:
        leaves[_occurrence_key(leaf.field, strip_derivatives=True)].append(leaf.node)

    placement = defaultdict(int)
    for field in completion.operator.indexed_fields:
        key = _occurrence_key(field, strip_derivatives=True)
        if not leaves[key]:
            raise UnsupportedAmplitudeAudit(
                f"Cannot match effective-operator field {field} to a partition leaf"
            )
        node = leaves[key].popleft()
        placement[node] += field.derivs
    return {node: degree for node, degree in placement.items() if degree}


def derivative_placement_polynomial(completion):
    """Return the formal external momentum polynomial for a completion.

    Each key is a sorted tuple ``((external_node, derivative_degree), ...)``.
    Momentum routed through a propagator cut is expanded multinomially over the
    leaves on one side.  An overall sign from choosing the opposite cut side is
    irrelevant to whether the complete tensor is zero.
    """

    contributions = tuple(getattr(completion, "momentum_contributions", ()))
    if not contributions:
        placement = _operator_derivative_placement(completion)
        return {tuple(sorted(placement.items())): 1}

    polynomial = {(): 1}
    for contribution in contributions:
        degree = contribution.derivative_degree
        if degree <= 0:
            continue
        side = contribution.cut_side or canonical_propagator_cut(
            completion.partition, completion.graph, contribution.edge
        )
        if not side:
            raise UnsupportedAmplitudeAudit(
                f"Empty propagator cut for edge {contribution.edge}"
            )
        factor = {}
        for allocation in weak_compositions(degree, len(side)):
            placement = tuple(
                (node, power)
                for node, power in zip(side, allocation)
                if power
            )
            factor[placement] = factor.get(placement, 0) + _multinomial(allocation)

        product = defaultdict(int)
        for left, left_coefficient in polynomial.items():
            for right, right_coefficient in factor.items():
                degrees = defaultdict(int, left)
                for node, power in right:
                    degrees[node] += power
                key = tuple(
                    sorted(
                        (node, power)
                        for node, power in degrees.items()
                        if power
                    )
                )
                product[key] += left_coefficient * right_coefficient
        polynomial = dict(product)

    expected_degree = sum(field.derivs for field in completion.operator.indexed_fields)
    actual_degrees = {sum(power for _, power in placement) for placement in polynomial}
    if actual_degrees != {expected_degree}:
        raise UnsupportedAmplitudeAudit(
            "Propagator momentum degree does not match the effective operator: "
            f"expected {expected_degree}, found {sorted(actual_degrees)}"
        )
    return polynomial


def audit_amplitude_symmetrisation(completion):
    """Return whether the fully induced formal gauge amplitude is identically zero."""

    try:
        polynomial = derivative_placement_polynomial(completion)
        total = 0
        witness = ""
        for placement, coefficient in sorted(polynomial.items()):
            amplitude = induced_amplitude_operator(
                completion, dict(placement)
            ).safe_simplify()
            if amplitude != 0 and not witness:
                witness = str(amplitude)
            total += coefficient * amplitude
        total = total.canon_bp() if total != 0 else 0
    except (AssertionError, KeyError, UnsupportedAmplitudeAudit, ValueError) as error:
        return AmplitudeAudit("unsupported", str(error), 0, "")

    if total == 0:
        return AmplitudeAudit(
            "zero",
            "Complete induced tensor vanishes after identical-field symmetrisation",
            len(polynomial),
            "",
        )
    return AmplitudeAudit("nonzero", "", len(polynomial), witness or str(total))
