#!/usr/bin/env python3

"""Explicit bases and exact projections for Lorentz-singlet contractions."""

from dataclasses import dataclass
from itertools import product
from math import gcd
from functools import reduce

from sympy import Matrix, Rational, ilcm

from neutrinomass.tensormethod.contract import lorentz_singlets
from neutrinomass.tensormethod.core import Index


LORENTZ_TYPES = ("u", "d")


def _perfect_pairings(ports):
    ports = tuple(ports)
    if not ports:
        return ((),)
    first = ports[0]
    pairings = []
    for position in range(1, len(ports)):
        second = ports[position]
        remaining = ports[1:position] + ports[position + 1 :]
        for tail in _perfect_pairings(remaining):
            pairings.append(((first, second), *tail))
    return tuple(pairings)


def _index_key(index):
    raised = index if index.is_up else -index
    return index.index_type, str(raised)


def _field_key(field):
    """Order external field slots without using Lorentz dummy-index names."""

    return (
        field.label,
        field.is_conj,
        field.derivs,
        field.dynkin,
        tuple(sorted((name, str(value)) for name, value in field.charges.items())),
        tuple((index.index_type, str(index)) for index in field.gauge_indices),
    )


def lorentz_field_key(field):
    """Return the stable external-field key used for Lorentz port ordering."""

    return _field_key(field)


def lorentz_field_port_layout(operator):
    """Return sorted external fields and their numbered Lorentz ports."""

    counts = {short_type: 0 for short_type in LORENTZ_TYPES}
    layout = []
    for field in sorted(operator.indexed_fields, key=_field_key):
        ports = {}
        for short_type in LORENTZ_TYPES:
            index_type = Index.get_index_types()[short_type]
            number = len(field.indices_by_type[index_type])
            start = counts[short_type]
            ports[short_type] = tuple(range(start, start + number))
            counts[short_type] += number
        layout.append((field, ports))
    return tuple(layout), tuple(
        (short_type, counts[short_type]) for short_type in LORENTZ_TYPES
    )


@dataclass(frozen=True)
class LorentzContraction:
    """A product of oriented two-spinor epsilon contractions."""

    coefficient: Rational
    pairings: tuple
    port_counts: tuple

    @property
    def label(self):
        structures = []
        for index_type, pairs in self.pairings:
            structures.append(
                index_type + ":" + ",".join(f"{left}-{right}" for left, right in pairs)
            )
        return "|".join(structures)

    def evaluation_vector(self):
        """Evaluate on all assignments of external spinors to a fixed 2D basis."""

        counts = dict(self.port_counts)
        offsets = {}
        total = 0
        for index_type in LORENTZ_TYPES:
            offsets[index_type] = total
            total += counts[index_type]

        values = []
        for assignment in product((0, 1), repeat=total):
            value = self.coefficient
            for index_type, pairs in self.pairings:
                offset = offsets[index_type]
                for left, right in pairs:
                    left_value = assignment[offset + left]
                    right_value = assignment[offset + right]
                    if left_value == right_value:
                        value = 0
                        break
                    value *= 1 if (left_value, right_value) == (0, 1) else -1
                if value == 0:
                    break
            values.append(value)
        return tuple(values)


def lorentz_contraction(operator, *, normalise=False):
    """Extract a dummy-name-independent epsilon-pairing representation."""

    fields = sorted(operator.indexed_fields, key=_field_key)
    index_to_port = {}
    port_counts = []
    for short_type in LORENTZ_TYPES:
        index_type = Index.get_index_types()[short_type]
        port = 0
        for field in fields:
            for index in field.indices_by_type[index_type]:
                key = _index_key(index)
                if key in index_to_port:
                    raise ValueError(
                        "A Lorentz index occurs on more than one field port"
                    )
                index_to_port[key] = (short_type, port)
                port += 1
        if port % 2:
            raise ValueError("A Lorentz singlet must have an even number of ports")
        port_counts.append((short_type, port))

    coefficient = Rational(operator.coeff)
    pairings = {index_type: [] for index_type in LORENTZ_TYPES}
    for invariant in operator.epsilons:
        first_type = invariant.indices[0].index_type
        short_type = next(
            (
                candidate
                for candidate in LORENTZ_TYPES
                if first_type == Index.get_index_types()[candidate]
            ),
            None,
        )
        if short_type is None:
            continue
        if len(invariant.indices) != 2:
            raise ValueError("Only two-index Lorentz invariants are supported")
        left_type, left = index_to_port[_index_key(invariant.indices[0])]
        right_type, right = index_to_port[_index_key(invariant.indices[1])]
        if left_type != short_type or right_type != short_type:
            raise ValueError("Lorentz invariant joins incompatible index types")
        if left > right:
            left, right = right, left
            coefficient *= -1
        pairings[short_type].append((left, right))

    for index_type, count in port_counts:
        if 2 * len(pairings[index_type]) != count:
            raise ValueError("Operator does not contain a complete Lorentz contraction")

    if normalise and coefficient:
        coefficient = Rational(1)
    return LorentzContraction(
        coefficient=coefficient,
        pairings=tuple(
            (index_type, tuple(sorted(pairings[index_type])))
            for index_type in LORENTZ_TYPES
        ),
        port_counts=tuple(port_counts),
    )


def primitive_coordinates(coordinates):
    """Remove an irrelevant common rational factor and fix the overall sign."""

    coordinates = tuple(Rational(value) for value in coordinates)
    nonzero = [value for value in coordinates if value]
    if not nonzero:
        return coordinates
    denominator = reduce(ilcm, (value.q for value in nonzero), 1)
    integers = [int(value * denominator) for value in coordinates]
    divisor = reduce(gcd, (abs(value) for value in integers if value))
    integers = [value // divisor for value in integers]
    first = next(value for value in integers if value)
    if first < 0:
        integers = [-value for value in integers]
    return tuple(Rational(value) for value in integers)


@dataclass(frozen=True)
class LorentzBasis:
    """A deterministic linearly independent Lorentz-singlet basis."""

    labels: tuple
    vectors: tuple
    port_counts: tuple

    @classmethod
    def from_operator(cls, operator):
        candidates = {}
        for singlet in lorentz_singlets(operator):
            if singlet.safe_simplify() == 0:
                continue
            contraction = lorentz_contraction(singlet, normalise=True)
            candidates[contraction.label] = contraction

        selected_labels = []
        selected_vectors = []
        rank = 0
        for label, contraction in sorted(candidates.items()):
            vector = contraction.evaluation_vector()
            trial = Matrix.hstack(
                *(Matrix(item) for item in (*selected_vectors, vector))
            )
            trial_rank = trial.rank()
            if trial_rank == rank:
                continue
            selected_labels.append(label)
            selected_vectors.append(vector)
            rank = trial_rank

        if not selected_vectors:
            raise ValueError("Operator has no non-zero Lorentz singlet")
        port_counts = next(iter(candidates.values())).port_counts
        return cls(tuple(selected_labels), tuple(selected_vectors), port_counts)

    @classmethod
    def from_port_counts(cls, port_counts):
        """Construct the full singlet basis without field-statistics reduction."""

        port_counts = tuple(port_counts)
        counts = dict(port_counts)
        candidates = {}
        pairing_options = [
            _perfect_pairings(range(counts[short_type]))
            for short_type in LORENTZ_TYPES
        ]
        for choices in product(*pairing_options):
            contraction = LorentzContraction(
                coefficient=Rational(1),
                pairings=tuple(zip(LORENTZ_TYPES, choices)),
                port_counts=port_counts,
            )
            candidates[contraction.label] = contraction

        selected_labels = []
        selected_vectors = []
        rank = 0
        for label, contraction in sorted(candidates.items()):
            vector = contraction.evaluation_vector()
            trial = Matrix.hstack(
                *(Matrix(item) for item in (*selected_vectors, vector))
            )
            trial_rank = trial.rank()
            if trial_rank == rank:
                continue
            selected_labels.append(label)
            selected_vectors.append(vector)
            rank = trial_rank
        if not selected_vectors:
            raise ValueError("Lorentz port content has no singlet")
        return cls(tuple(selected_labels), tuple(selected_vectors), port_counts)

    @property
    def dimension(self):
        return len(self.labels)

    def project(self, operator, *, primitive=True):
        contraction = lorentz_contraction(operator)
        if contraction.port_counts != self.port_counts:
            raise ValueError("Operator has different Lorentz field content")
        basis_matrix = Matrix.hstack(*(Matrix(vector) for vector in self.vectors))
        target = Matrix(contraction.evaluation_vector())
        try:
            coordinates, parameters = basis_matrix.gauss_jordan_solve(target)
        except ValueError:
            return None
        if parameters.rows:
            raise ValueError("Lorentz basis is not linearly independent")
        result = tuple(Rational(value) for value in coordinates)
        return primitive_coordinates(result) if primitive else result
