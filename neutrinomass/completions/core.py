#!/usr/bin/env python3

"""Core classes and functions for completions code."""

from typing import Dict
from copy import deepcopy

from functools import reduce

from neutrinomass.tensormethod import IndexedField
from neutrinomass.tensormethod.core import BOSE, FERMI, Index, eps, get_dynkin

# from neutrinomass.completions.tikzfeynman import tikz_export


class FieldType(IndexedField):
    """Base class for exotic fields."""

    def __new__(cls, *args, **kwargs):
        return super(FieldType, cls).__new__(cls, *args, **kwargs)

    def lower_su2(self):
        undotted, dotted, _, isospin, _ = self.indices_by_type.values()
        epsilons = []
        partner = self
        for idx in [*undotted, *dotted, *isospin]:
            lower = str(idx) + "^"
            partner = partner.substituted_indices((idx, lower))
            epsilon = eps("-" + lower + " -" + str(idx))
            epsilons.append(epsilon)

        return cons_completion_field(partner), epsilons

    @property
    def indexed_field(self):
        return IndexedField(label=self.label, indices=self.index_labels)

    # @property
    # def info(self):
    #     return {
    #         "label": self.label,
    #         "indices": self.indices,
    #         "charges": self.charges,
    #         "is_conj": self.is_conj,
    #         "symmetry": self.symmetry,
    #         "latex": self.latex,
    #     }

    def __hash__(self):
        dict_ = {
            "indices": tuple(sorted(self.indices)),
            "label": self.label,
            "dynkin": self.dynkin,
            "charges": tuple(sorted(self.charges.items())),
        }
        return hash(tuple(dict_.items()))

    def __deepcopy__(self, memo):
        return self.__class__(
            label=self.label,
            indices=deepcopy(self.indices, memo),
            charges=deepcopy(self.charges, memo),
            latex=self.latex,
            is_conj=self.is_conj,
            symmetry=deepcopy(self.symmetry, memo),
            comm=self.comm,
        )


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
            # **kwargs,
        )

        assert self.is_scalar


def assert_real_rep(indices: str, charges) -> None:
    if charges and "y" in charges:
        assert charges["y"] == 0

    colour_dynkin_str = get_dynkin(indices)[2:4]

    # simplistic way to check if colour rep is real
    assert colour_dynkin_str == "".join(reversed(colour_dynkin_str))


class RealScalar(FieldType):
    def __init__(
        self, label, indices, charges=None, latex=None, is_conj=False, **kwargs
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
            # **kwargs,
        )

        assert self.is_scalar
        assert_real_rep(indices, charges)


class MajoranaFermion(FieldType):
    def __init__(
        self, label, indices, charges=None, latex=None, is_conj=False, **kwargs
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
            comm=FERMI,
            latex=latex,
            # **kwargs,
        )

        assert self.is_fermion
        assert_real_rep(indices, charges)

    def majorana_partner(self):
        """New copy of fermion with colour indices flipped"""
        undotted, dotted, colour, isospin, _ = Index.indices_by_type(
            self.indices
        ).values()
        colour = tuple(i.conj for i in colour)
        indices = undotted + dotted + colour + isospin

        return MajoranaFermion(
            self.label,
            indices,
            latex=self.latex,
            is_conj=self.is_conj,
            symmetry=self.symmetry,
            comm=FERMI,
            charges=None,
        )


class VectorLikeDiracFermion(FieldType):
    """Stands in for two fermion fields distinguished here only by dotted and
    undotted indices.

    Example:
        >>> psi = VectorLikeDiracFermion("ψ", "u0 i1")
        >>> psi
        ψ(u0, i1)
        >>> psi.dirac_partner()
        ψ~(u0, i1)

    Here ψ(u0, i1) and ψ~(u0, i1) are different fields.

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
            latex=latex,
            # **kwargs,
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
        symb = "~"  # indicate bar with circumflex
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

    def __hash__(self):
        return hash((self.name, self.operator.simplify()))


class FailedCompletion:
    def __init__(self, reason: str = ""):
        self.reason = reason


class Completion:
    def __init__(self, operator, partition, graph, exotics, terms):
        self.operator = operator
        self.partition = partition
        self.graph = graph
        self.exotics = exotics
        self.terms = terms

    def __eq__(self, other):
        if not isinstance(other, Completion):
            return False
        return (
            self.operator == other.operator
            and self.exotic_info() == other.exotic_info()
            and self.partition == other.partition
        )

    def __hash__(self):
        return hash((self.operator, self.exotic_info(), self.partition))

    @property
    def lagrangian(self):
        return Lagrangian(terms=self.terms)

    def draw_diagram(self):
        import matplotlib.pyplot as plt
        import networkx as nx

        g = self.graph
        edge_labels = nx.get_edge_attributes(g, name="particle")
        pos = nx.spring_layout(g)

        plt.figure()
        nx.draw(g, pos=pos, edge_color="black", node_size=0)
        nx.draw_networkx_edge_labels(
            g, pos=pos, edge_labels=edge_labels, font_color="red"
        )
        plt.axis("off")
        plt.show()

    def exotic_info(self) -> Dict[FieldType, tuple]:
        info = {}
        for e in self.exotics:
            # normalise hypercharge to be positive
            if e.y < 0:
                charges = sorted(e.conj.charges.items())
                sm = e.conj.sm_irrep
            else:
                charges = sorted(e.charges.items())
                sm = e.sm_irrep

            if e.is_fermion:
                lorentz = "F"
            elif e.is_scalar:
                lorentz = "S"
            else:
                raise ValueError("Unrecognised exotic field type.")

            info[e] = (lorentz,) + sm + tuple(charges)

        return info

    def exotic_fields(self):
        return set([e.field for e in self.exotics])

    def export(self, path: str):
        # currently gives you circular import

        # export graph tex file
        # tikz_export(g=self.graph, path=path)
        pass


def cons_completion_field(indexed_field: IndexedField) -> FieldType:
    label = indexed_field.label
    indices = indexed_field.indices
    charges = indexed_field.charges
    latex = indexed_field.latex

    if indexed_field.is_fermion:
        if indexed_field.is_real_sm_irrep:
            return MajoranaFermion(label, indices, charges=charges, latex=latex)

        return VectorLikeDiracFermion(label, indices, charges=charges, latex=latex)

    if indexed_field.is_scalar:
        if indexed_field.is_real_sm_irrep:
            return RealScalar(label, indices, charges=charges, latex=latex)

        return ComplexScalar(label, indices, charges=charges, latex=latex)

    raise Exception("Unrecognised Lorentz structure in field.")
