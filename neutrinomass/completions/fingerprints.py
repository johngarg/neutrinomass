#!/usr/bin/env python3

"""Deterministic fingerprints for completion-regression baselines."""

from hashlib import sha256
from typing import Iterable, Tuple

import networkx as nx

from neutrinomass.completions.core import Completion
from neutrinomass.completions.equivalence import (
    conjugated_interaction_graph,
    interaction_graph,
)


def interaction_fingerprint(term) -> str:
    """Return a stable structural hash of an interaction contraction graph.

    Exact equivalence decisions use graph isomorphism in ``equivalence.py``;
    this collision-resistant digest is only a compact regression identifier.
    """

    def graph_hash(graph):
        graph = graph.copy()
        nx.set_node_attributes(
            graph,
            {
                node: sha256(
                    repr(data["signature"]).encode("utf-8")
                ).hexdigest()
                for node, data in graph.nodes(data=True)
            },
            "colour",
        )
        return nx.weisfeiler_lehman_graph_hash(
            graph,
            node_attr="colour",
            iterations=max(3, graph.number_of_nodes()),
            digest_size=32,
        )

    return min(
        graph_hash(interaction_graph(term)),
        graph_hash(conjugated_interaction_graph(term)),
    )


def democratic_model_fingerprint(completion: Completion) -> tuple:
    """Return the unique particle-species content used by democratic filtering."""

    return tuple(sorted(set(completion.exotic_info().values())))


def lagrangian_fingerprint(completion: Completion) -> Tuple[str, ...]:
    """Return sorted structural interaction hashes for regression snapshots."""

    return tuple(sorted(interaction_fingerprint(term) for term in completion.terms))


def completion_fingerprint(completion: Completion) -> tuple:
    """Return the deterministic physics and provenance fingerprint of a completion."""

    derivative_edges = tuple(
        sorted(tuple(sorted(edge)) for edge in completion.derivative_edges)
    )
    projection = getattr(completion, "lorentz_projection", None)
    projection_fingerprint = None
    if projection is not None:
        projection_fingerprint = (
            projection.basis_labels,
            projection.coordinates,
            projection.derivative_field,
        )
    fingerprint = (
        completion.operator.name,
        completion.topology,
        completion.canonical_topology,
        democratic_model_fingerprint(completion),
        lagrangian_fingerprint(completion),
        derivative_edges,
    )
    if projection_fingerprint is not None:
        return fingerprint + (projection_fingerprint,)
    return fingerprint


def completion_digest(completions: Iterable[Completion]) -> str:
    """Hash a completion multiset without depending on generation order."""

    return completion_fingerprint_digest(
        completion_fingerprint(completion) for completion in completions
    )


def completion_fingerprint_digest(fingerprints: Iterable[tuple]) -> str:
    """Hash precomputed completion fingerprints without retaining completions."""

    payload = "\n".join(sorted(repr(fingerprint) for fingerprint in fingerprints))
    return sha256(payload.encode()).hexdigest()
