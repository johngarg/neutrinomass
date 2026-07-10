#!/usr/bin/env python3

"""Deterministic fingerprints for completion-regression baselines."""

from hashlib import sha256
from typing import Iterable, Tuple

from neutrinomass.completions.core import Completion
from neutrinomass.tensormethod.utils import safe_nocoeff


def interaction_fingerprint(term) -> str:
    """Return an exact, coefficient-free tensor fingerprint for one interaction."""

    canonical_term = term.fill_free_indices().safe_simplify()
    return str(safe_nocoeff(canonical_term))


def democratic_model_fingerprint(completion: Completion) -> tuple:
    """Return the unique particle-species content used by democratic filtering."""

    return tuple(sorted(set(completion.exotic_info().values())))


def lagrangian_fingerprint(completion: Completion) -> Tuple[str, ...]:
    """Return the sorted exact interaction fingerprints of a completion."""

    return tuple(sorted(interaction_fingerprint(term) for term in completion.terms))


def completion_fingerprint(completion: Completion) -> tuple:
    """Return the deterministic physics and provenance fingerprint of a completion."""

    derivative_edges = tuple(
        sorted(tuple(sorted(edge)) for edge in completion.derivative_edges)
    )
    return (
        completion.operator.name,
        completion.topology,
        completion.canonical_topology,
        democratic_model_fingerprint(completion),
        lagrangian_fingerprint(completion),
        derivative_edges,
    )


def completion_digest(completions: Iterable[Completion]) -> str:
    """Hash a completion multiset without depending on generation order."""

    payload = "\n".join(
        sorted(repr(completion_fingerprint(completion)) for completion in completions)
    )
    return sha256(payload.encode()).hexdigest()
