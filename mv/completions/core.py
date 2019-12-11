#!/usr/bin/env python3

"""Core classes and functions for completions code."""

from functools import reduce

from mv.tensormethod import IndexedField
from mv.tensormethod.core import BOSE, FERMI, Index


class FieldType(IndexedField):
    """Base class for exotic fields."""

    def __new__(cls, *args, **kwargs):
        return super(FieldType, cls).__new__(cls, *args, **kwargs)

    def __mul__(self, other):
        if isinstance(other, self.__class__):
            return TensorProduct(self, other)
        elif isinstance(other, TensorProduct):
            return TensorProduct(self, *other.tensors)


class ComplexScalar(FieldType):
    def __init__(
        self,
        label,
        indices,
        charges=None,
        latex=None,
        is_conj=False,
        symmetry=None,
        **kwargs,
    ):
        if isinstance(indices, list):
            indices = " ".join(str(i) for i in indices)

        IndexedField.__init__(
            self,
            label=label,
            indices=indices,
            charges=charges,
            is_conj=is_conj,
            symmetry=None,
            comm=BOSE,
            latex=latex,
        )

        assert self.is_scalar


class RealScalar(FieldType):
    def __init__(self, label, indices, latex=None):
        if isinstance(indices, list):
            indices = " ".join(str(i) for i in indices)

        IndexedField.__init__(
            self,
            label=label,
            indices=indices,
            charges=None,
            is_conj=False,
            symmetry=None,
            comm=BOSE,
            latex=latex,
        )

        assert self.is_scalar


class MajoranaFermion(FieldType):
    def __init__(self, label, indices, latex=None):
        if isinstance(indices, list):
            indices = " ".join(str(i) for i in indices)

        IndexedField.__init__(
            self,
            label=label,
            indices=indices,
            charges=None,
            is_conj=False,
            symmetry=None,
            comm=FERMI,
            latex=latex,
        )

        assert self.is_fermion


class VectorLikeDiracFermion(FieldType):
    """Stands in for two fermion fields distinguished here only by dotted and
    undotted indices.

    Example:
        >>> psi = VectorLikeDiracFermion("ψ", "u0 i1")
        >>> psi
        ψ(u0, i1)
        >>> psi.dirac_partner()
        ψ(d0, i1)

    Here ψ(u0, i1) and ψ(d0, i1) are different fields.

    """

    def __init__(
        self,
        label,
        indices,
        charges=None,
        latex=None,
        is_unbarred=True,
        is_conj=False,
        symmetry=None,
        comm=FERMI,
    ):
        if isinstance(indices, list):
            indices = " ".join(str(i) for i in indices)

        IndexedField.__init__(
            self,
            label=label,
            indices=indices,
            charges=charges,
            is_conj=is_conj,
            symmetry=None,
            latex=latex,
            **kwargs,
        )

        assert self.is_fermion
        self.is_unbarred = is_unbarred

    @property
    def conj(self):
        """Returns a copy of self but conjugated"""

        is_conj = self.is_conj
        if is_conj:
            label = self.label.replace("†", "")
        else:
            label = self.label + "†"

        return self.__class__(
            label=label,
            indices=" ".join(i.conj.label for i in self.indices),
            charges={k: -v for k, v in self.charges.items()},
            is_conj=(not is_conj),
            symmetry=self.symmetry,
            comm=self.comm,
            is_unbarred=self.is_unbarred,
        )

    def dirac_partner(self):
        undotted, dotted, colour, isospin, _ = Index.indices_by_type(
            self.indices
        ).values()
        colour = tuple(i.conj for i in colour)
        indices = undotted + dotted + colour + isospin
        charges = {k: -v for k, v in self.charges.items()}

        is_unbarred = self.is_unbarred
        symb = "^"  # indicate bar with circumflex
        if is_unbarred:
            label = self.label + symb
        else:
            label = self.label.replace(symb, "")

        return VectorLikeDiracFermion(
            label,
            indices,
            charges=charges,
            latex=self.latex,
            is_unbarred=(not self.is_unbarred),
            is_conj=self.is_conj,
            symmetry=self.symmetry,
            comm=FERMI,
        )


class EffectiveOperator:
    def __init__(self, name, operator):
        self.name = name
        self.operator = operator

    @property
    def fields(self):
        return self.operator.fields

    @property
    def indexed_fields(self):
        return self.operator.indexed_fields

    @property
    def mass_dimension(self):
        d = sum(f.mass_dim for f in self.fields)
        d_int = int(d)
        assert d == d_int
        return d_int

    @property
    def topology_type(self):
        """Returns a dictionary {"n_scalars": n_scalars, "n_fermions": n_fermions}."""

        n_scalars, n_fermions = 0, 0
        for f in self.fields:
            if f.is_scalar:
                n_scalars += 1
            elif f.is_fermion:
                n_fermions += 1

        return {"n_scalars": n_scalars, "n_fermions": n_fermions}


class Completion:
    def __init__(self, partition, graph, exotics):
        self.partition = partition
        self.graph = graph

    @property
    def lagrangian(self):
        pass

    @property
    def diagram(self):
        pass


def cons_completion_field(indexed_field: IndexedField) -> FieldType:
    label = indexed_field.label
    indices = indexed_field.indices
    charges = indexed_field.charges
    latex = indexed_field.latex

    if indexed_field.is_fermion:
        if indexed_field.is_real_sm_irrep:
            return MajoranaFermion(label, indices, latex=latex)

        return VectorLikeDiracFermion(label, indices, charges=charges, latex=latex)

    if indexed_field.is_scalar:
        if indexed_field.is_real_sm_irrep:
            return RealScalar(label, indices, latex=latex)

        return ComplexScalar(label, indices, charges=charges, latex=latex)

    raise Exception("Unrecognised Lorentz structure in field.")
