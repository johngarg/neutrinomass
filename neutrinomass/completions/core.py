#!/usr/bin/env python3

"""Core classes and functions for completions code."""

import sys
from typing import Dict
from copy import deepcopy

from functools import reduce

from neutrinomass.tensormethod.core import (
    BOSE,
    FERMI,
    Index,
    IndexedField,
    eps,
    get_dynkin,
)
from neutrinomass.tensormethod.lagrangian import Lagrangian

# from neutrinomass.completions.tikzfeynman import tikz_export


class FieldType(IndexedField):
    """Base class for exotic fields."""

    def __new__(cls, *args, **kwargs):
        return super(FieldType, cls).__new__(cls, *args, **kwargs)

    def lower_su2(self, skip=[]):
        undotted, dotted, _, isospin, _ = self.indices_by_type.values()
        epsilons = []
        partner = self
        for idx in [*undotted, *dotted, *isospin]:
            if idx.index_type in skip:
                continue
            lower = str(idx) + "^"
            partner = partner.substituted_indices((idx, lower))
            epsilon = eps("-" + lower + " -" + str(idx))
            epsilons.append(epsilon)

        # For VectorLikeDiracFermion, keep track of (un)barred boolean
        is_unbarred = None
        if hasattr(self, "is_unbarred"):
            is_unbarred = self.is_unbarred

        return cons_completion_field(partner, is_unbarred=is_unbarred), epsilons

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

    @property
    def mass_term(self):
        lower, epsilons = self.conj.lower_su2()
        return reduce(lambda x, y: x * y, epsilons, lower * self)

    # @property
    # def kinetic_term(self):
    #     """I think this will mess with U(1) symmetries algorithm as currently
    #     implemented"""
    #     lower, epsilons = self.conj.lower_su2()
    #     all_og_indices = "u0 d0 " + " ".join(str(i) for i in self.indices)
    #     all_lower_indices = "u1 d1 " + " ".join(str(i) for i in lower.indices)
    #     deriv_og = D(lower, "11")(all_og_indices)
    #     deriv_lower = D(lower, "11")(all_lower_indices)
    #     term = deriv_og * deriv_lower
    #     return reduce(lambda x, y: x * y, epsilons, term)


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

    def swap_colour_indices(self):
        """New copy of field with colour indices flipped.

        # TODO Refactor this out with majorana_partner to FieldType
        """
        undotted, dotted, colour, isospin, _ = Index.indices_by_type(
            self.indices
        ).values()
        colour = tuple(i.conj for i in colour)
        indices = undotted + dotted + colour + isospin

        return RealScalar(
            self.label,
            indices,
            latex=self.latex,
            is_conj=self.is_conj,
            symmetry=self.symmetry,
            comm=BOSE,
            charges=None,
        )

    @property
    def mass_term(self):
        lower, epsilons = self.swap_colour_indices().lower_su2()
        return reduce(lambda x, y: x * y, epsilons, lower * self)


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

    @property
    def conj_indices(self):
        return self.conj

    @property
    def mass_term(self):
        lower, epsilons = self.majorana_partner().lower_su2()
        return reduce(lambda x, y: x * y, epsilons, lower * self)


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
            comm=FERMI,
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
            is_unbarred=self.is_unbarred,
            symmetry=self.symmetry,
            comm=self.comm,
        )

    @property
    def conj_indices(self):
        """Returns a copy of self conjugated but leaves charges alone."""
        is_conj = self.is_conj
        if is_conj:
            label = self.label.replace("†", "")
        else:
            label = self.label + "†"

        return self.__class__(
            label=label,
            indices=" ".join(i.conj.label for i in self.indices),
            charges=self.charges,
            is_conj=(not is_conj),
            is_unbarred=self.is_unbarred,
            symmetry=self.symmetry,
            comm=self.comm,
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

    @property
    def mass_term(self):
        lower, epsilons = self.dirac_partner().lower_su2()
        return reduce(lambda x, y: x * y, epsilons, lower * self)


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
        d += sum(f.derivs for f in self.fields)
        d_int = int(d)
        assert d == d_int
        return d_int

    @property
    def topology_type(self):
        """Returns a dictionary {"n_scalars": n_scalars, "n_fermions": n_fermions}."""

        n_scalars, n_fermions = 0, 0
        for f in self.fields:
            if f.is_boson:
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
    def __init__(self, operator, partition, graph, exotics, terms, topology=None):
        self.operator = operator
        self.partition = partition
        self.graph = graph
        self.exotics = exotics
        self.terms = terms
        self.topology = topology

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

    def __deepcopy__(self, memo):
        return self.__class__(
            operator=deepcopy(self.operator, memo),
            partition=deepcopy(self.partition, memo),
            graph=deepcopy(self.graph, memo),
            exotics=deepcopy(self.exotics, memo),
            terms=deepcopy(self.terms, memo),
        )

    @property
    def lagrangian(self):
        return Lagrangian(exotics=self.exotics, interaction_terms=self.terms)

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

    def info(self):
        print("Fields:")
        for k, v in self.exotic_info().items():
            print("{:<5s}{:<20s}".format(k.label, format_quantum_numbers(v)))

        print("\nLagrangian:")
        in_jupyter = sys.argv[-1].endswith("json")
        for term in self.terms:
            if in_jupyter:
                display(term)
            else:
                print(term)

        print("\nDiagram:")
        self.draw_diagram()


class Model:
    def __init__(self, completions):
        self.completions = completions

    @property
    def exotic_numbers(self):
        return sorted(
            set(
                format_quantum_numbers(i)
                for i in self.completions[0].exotic_info().values()
            )
        )

    def __repr__(self):
        return "Model(" + " + ".join(i for i in self.exotic_numbers) + ")"


def format_quantum_numbers(info: tuple):
    """Takes an expression like

    ('S', 1, 0, 2, ('3b', 1), ('y', 2/3))

    and returns a string like

    S(3, 3, 2/3)(3b: 1)

    """
    lorentz, su3_up, su3_down, su2, *charges = info
    su3_dim = lambda m, n: 0.5 * (m + 1) * (n + 1) * (m + n + 2)
    # For now just add the bar for more lowered than raised indices, but for
    # larger reps this will be problematic
    su3_dim_format = lambda m, n: str(int(su3_dim(m, n))) + ("b" if n > m else "")
    charges_dict = dict(charges)
    return f"{lorentz}({su3_dim_format(int(su3_up), int(su3_down))}, {str(int(su2) + 1)}, {charges_dict['y']})({charges_dict['3b']})"


def cons_completion_field(indexed_field: IndexedField, is_unbarred=None) -> FieldType:
    label = indexed_field.label
    indices = indexed_field.indices
    charges = indexed_field.charges
    latex = indexed_field.latex
    is_conj = indexed_field.is_conj

    if indexed_field.is_fermion:
        if indexed_field.is_real_sm_irrep:
            return MajoranaFermion(
                label, indices, charges=charges, latex=latex, is_conj=is_conj
            )

        if is_unbarred is None:
            is_unbarred = True

        return VectorLikeDiracFermion(
            label,
            indices,
            charges=charges,
            latex=latex,
            is_unbarred=is_unbarred,
            is_conj=is_conj,
        )

    if indexed_field.is_scalar:
        if indexed_field.is_real_sm_irrep:
            return RealScalar(label, indices, charges=charges, latex=latex)

        return ComplexScalar(
            label, indices, charges=charges, latex=latex, is_conj=is_conj
        )

    raise Exception("Unrecognised Lorentz structure in field.")
