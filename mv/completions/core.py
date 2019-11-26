#!/usr/bin/env python3

"""Core classes and functions for completions code."""

from functools import reduce

from mv.tensormethod import IndexedField
from mv.tensormethod.core import BOSE, FERMI


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
    def __init__(self, label, indices, charges, latex=None):
        IndexedField.__init__(
            self,
            label=label,
            indices=indices,
            charges=charges,
            is_conj=False,
            symmetry=None,
            comm=BOSE,
            latex=latex,
        )

        assert self.is_scalar


class RealScalar(FieldType):
    def __init__(self, label, indices, latex=None):
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
    def __init__(self, label, indices, charges=None, latex=None):
        IndexedField.__init__(
            self,
            label=label,
            indices=indices,
            charges=charges,
            is_conj=False,
            symmetry=None,
            comm=FERMI,
            latex=latex,
        )

        assert self.is_fermion


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
    def __init__(self, partition, graph):
        self.partition = partition
        self.graph = graph

    @property
    def lagrangian(self):
        pass

    @property
    def diagram(self):
        pass
