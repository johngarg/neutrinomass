#!/usr/bin/env python


from collections import defaultdict
from copy import copy, deepcopy
from itertools import groupby
from string import ascii_lowercase
from typing import Iterable
from typing import List
from typing import NamedTuple
from typing import Tuple
from typing import Union
from functools import reduce

import sympy.tensor.tensor as tensor
from basisgen import irrep
from sympy import flatten, Rational
from sympy.core.numbers import Zero

from neutrinomass.tensormethod.lnv import BL_LIST
from neutrinomass.tensormethod.utils import (
    repr_tree,
    to_tex,
    TEX_GREEK_LOWERCASE,
    DOTTED_TEX_GREEK_LOWERCASE,
    strip_parens,
)

FERMI, BOSE = "fermi", "bose"
tensor.TensorManager.set_comms(
    (FERMI, FERMI, 1), (BOSE, BOSE, 0), (FERMI, BOSE, 0), (BOSE, FERMI, 0)
)

CHARGES = ("y", "3b")


class History(NamedTuple):
    """A structure representing the left and right parents of the result of a binary
    operation.

    """

    left: "Field"
    right: "Field"


class Prod(NamedTuple):
    """A structure representing the symbolic product of two fields.

    Used when multiplying fields to keep track of multiplications.

    Example:
        >>> Prod(irrep=LL(00000), left=L(10001), right=L(10001))

    """

    irrep: "Field"
    left: Union["Prod", "Field"]
    right: Union["Prod", "Field"]


# Use sympy tensor indices as base for custom index types
ISOSPIN = tensor.TensorIndexType("Isospin", metric=True, dummy_fmt="I", dim=2)
COLOUR = tensor.TensorIndexType("Colour", metric=None, dummy_fmt="C", dim=3)
GENERATION = tensor.TensorIndexType("Generation", metric=None, dummy_fmt="G", dim=3)
UNDOTTED = tensor.TensorIndexType("Undotted", metric=True, dummy_fmt="U", dim=2)
DOTTED = tensor.TensorIndexType("Dotted", metric=True, dummy_fmt="D", dim=2)


class Index(tensor.TensorIndex):
    """A tensor index."""

    def __new__(cls, label: str, *args, **kwargs):
        if isinstance(label, Index):
            return label

        # This will be called in __mul__ method of IndexedTensor
        # i.e. Will be a contracted index
        # if isinstance(label, tensor.TensorIndex):
        label = str(label).lower()

        # deal with negative sign for down index
        is_up = True
        if label[0] == "-":
            is_up = False
            label = label[1:]

        dummy = True if label[0] == "_" else False
        if dummy:
            label = label[1:] + "_"

        tensor_type = cls.classify_index(label)
        return super(Index, cls).__new__(
            cls, name=label, tensortype=tensor_type, is_up=is_up
        )

    def __init__(self, label, *args, **kwargs):
        if isinstance(label, Index):
            self.label = label.label
        self.label = str(label)  # the label that was passed in

    def __neg__(self):
        is_up = self.is_up
        new_label = "-" + self.label if is_up else self.label[1:]
        return Index(label=new_label)

    @property
    def conj(self):
        # conj takes you between dotted and undotted
        if self.index_type in Index.get_lorentz_index_types().values():
            if self.index_type == "Undotted":
                return Index(label=self.label.replace("u", "d"))
            return Index(label=self.label.replace("d", "u"))

        # conj does nothing to isospin indices (implicit epsilon)
        # conj also does nothing to generation indices
        if self.index_type == "Isospin" or self.index_type == "Generation":
            return self

        # generally lower index (only affects colour here)
        return -self

    @property
    def index_type(self):
        return str(self.tensor_index_type)

    @classmethod
    def dynkin(cls, indices):
        """Returns a tuple (# raised, # lowered) indices.

        indices: collection of Index objects

       """
        up, down = 0, 0
        for i in indices:
            if i.is_up:
                up += 1
            else:
                down += 1
        return up, down

    @classmethod
    def get_dynkin_labels(cls):
        """[(is_raised, type), ...]"""
        return [(True, "u"), (True, "d"), (True, "c"), (False, "c"), (True, "i")]

    @classmethod
    def get_tensor_index_types(cls):
        """Map between shorthand index-type labels and the sympy tensor index types."""
        return {"u": UNDOTTED, "d": DOTTED, "c": COLOUR, "i": ISOSPIN, "g": GENERATION}

    @classmethod
    def get_index_types(cls):
        """Map between shorthand index-type labels and names of index types."""
        return {
            "u": "Undotted",
            "d": "Dotted",
            "c": "Colour",
            "i": "Isospin",
            "g": "Generation",
        }

    @classmethod
    def get_index_labels(cls):
        """Map between names of index types and shorthand index-type labels."""
        return {v: k for k, v in cls.get_index_types().items()}

    @classmethod
    def fresh(cls, type_) -> "Index":
        """Return a fresh, unused index by mutating a global counter internal to sympy's
        tensor module.

        """
        idx = cls(tensor.TensorIndex(True, cls.get_tensor_index_types()[type_]))
        label = idx.label.replace("i", type_)
        return cls(label)

    @classmethod
    def fresh_indices(cls, dynkin_str) -> str:
        """Returns a string of fresh index names matching ``dynkin_str``"""
        out = []
        for (raised, type_), number in zip(cls.get_dynkin_labels(), dynkin_str):
            for _ in range(int(number)):
                idx = Index.fresh(type_) if raised else -Index.fresh(type_)
                out.append(idx)
        return " ".join(i.__repr__() for i in out)

    @classmethod
    def get_sm_index_types(cls):
        return {"c": "Colour", "i": "Isospin"}

    @classmethod
    def get_lorentz_index_types(cls):
        return {"u": "Undotted", "d": "Dotted"}

    @property
    def is_su2(self):
        type_ = self.index_type
        return type_ == "Isospin" or type_ == "Undotted" or type_ == "Dotted"

    @classmethod
    def classify_index(cls, idx: str) -> tensor.TensorIndexType:
        char = idx[0] if idx[0] != "-" else idx[1]
        return Index.get_tensor_index_types()[char]

    @classmethod
    def indices_by_type(cls, indices):
        """Returns a dictionary mapping index type to tuple of indices of that type."""
        result = {k: [] for k in Index.get_index_types().values()}
        for i in indices:
            result[i.index_type].append(i)
        return {k: tuple(v) for k, v in result.items()}

    @classmethod
    def cons_index(cls, index_type: str, label: int, is_up: bool) -> "Index":
        """Construct an index with provided properties."""
        prefix = Index.get_index_labels()[index_type]
        prefix = ("-" if not is_up else "") + prefix
        return Index(prefix + str(label))

    @classmethod
    def conj_index_string(cls, index_string: str) -> str:
        """Return string of conjugated indices"""
        indices = [cls(s).conj for s in index_string.split(" ")]
        return " ".join(str(i) for i in indices)

    @property
    def _dict(self):
        return {"label": self.label}


class Field:
    def __init__(
        self,
        label: str,  # the str representation of the field
        dynkin: str,  # a str with the dynkin digits
        charges=None,  # a dictionary mapping charges to real values
        is_conj=False,  # conjugated flag
        comm=0,  # commutation symmetry
        symmetry=None,  # tensor index symmetry
        history=None,  # parents
        multiplicity=1,  # number of possibilities
        latex=None,  # base latex representation
        nf=1,  # number of generations (currently 1 or 3)
        derivs=0,  # number of derivatives acting on field
        stripped=None,  # field stripped of a derivative
        **kwargs,
    ):
        """A representation of a tensor transforming under Lorentz x SM with charges and
        a label.

        Tensor symmetry inherited from sympy, possible options:
            ``[[1]]``         vector
            ``[[1]*n]``       symmetric tensor of rank ``n``
            ``[[n]]``         antisymmetric tensor of rank ``n``
            ``[[2, 2]]``      monoterm slot symmetry of the Riemann tensor
            ``[[1],[1]]``     vector*vector
            ``[[2],[1],[1]]`` (antisymmetric tensor)*vector*vector

        Example:
            >>> Field("A", dynkin="10011", charges={"y": 1})
            A(10011)(1)

        """

        # initialise charges
        default_charges = {k: 0 for k in CHARGES}  # {"y": 0, "3b": 0}
        if charges is None:
            charges = default_charges
        else:
            # Will fill in absent keys with default values
            charges = {**default_charges, **charges}

        # make sure charges contains hypercharge
        assert "y" in charges.keys()

        # allow charges to be passed in as strings
        proc_charges = {}
        for k, v in charges.items():
            proc_charges[k] = Rational(v)
        charges = proc_charges

        # for SU(2) and SU(3), indices will always be symmetric
        if symmetry is None:
            symmetry = []
            for i in map(int, dynkin):
                if i:
                    symmetry += [[1] * i]

        if not isinstance(dynkin, str):
            dynkin = "".join(str(i) for i in dynkin)

        if latex is None:
            latex = to_tex(label)

        if nf != 1 and nf != 3:
            raise ValueError("Only currently supporting 1 or 3 generations.")

        self.symmetry = symmetry
        self.history = history
        self.charges = charges
        self.multiplicity = multiplicity
        self.label = label
        self.dynkin = dynkin.strip()
        self.comm = comm
        self.is_conj = is_conj
        self.latex = latex
        self.nf = nf
        self.derivs = derivs
        self.stripped = stripped

    def __call__(self, indices: str) -> "IndexedField":
        """Returns an IndexedField object.

        Indices must match dynkin structure.

        """
        dynkin_ints = [int(i) for i in self.dynkin]
        su2_plus, su2_minus, su3_up, su3_down, su2 = dynkin_ints
        index_types = (
            [UNDOTTED] * su2_plus
            + [DOTTED] * su2_minus
            + [COLOUR] * (su3_up + su3_down)
            + [ISOSPIN] * su2
        )

        # make sure names of indices are correct (to avoid user confusion)
        assert_consistent_indices(
            [i for i in indices.split() if i[0] != "g"],  # ignore generation indices
            index_types,
            (su3_up, su3_down),
        )
        return IndexedField(
            label=self.label_with_dagger,
            indices=indices,
            charges=self.charges,
            is_conj=self.is_conj,
            symmetry=self.symmetry,
            # multiplicity=self.multiplicity,
            comm=self.comm,
            latex=self.latex,
            nf=self.nf,
            derivs=self.derivs,
            stripped=self.stripped,
        )

    @property
    def label_with_dagger(self):
        maybe_conj = "†" if self.is_conj else ""
        return self.label + maybe_conj

    def __repr__(self):
        return self.label_with_dagger + f"({self.dynkin})({self.y})"

    __str__ = __repr__

    @property
    def qn_string(self):
        if self.is_scalar:
            lor = "S"
        elif self.is_fermion:
            lor = "F"
        else:
            lor = "X"

        return lor + self.dynkin[2:] + f"({self.y})({self.charges['3b']})"

    @property
    def y(self):
        """Return the hypercharge as a sympy rational object"""
        return self.charges["y"]

    @property
    def dynkin_ints(self):
        return tuple([int(i) for i in self.dynkin])

    @property
    def lorentz_irrep(self):
        return self.dynkin_ints[:2]

    @property
    def sm_irrep(self):
        return self.dynkin_ints[2:]

    @property
    def colour_irrep(self):
        return self.sm_irrep[:2]

    @property
    def isospin_irrep(self):
        return self.sm_irrep[2:]

    @property
    def quantum_numbers(self):
        return self.sm_irrep + (self.y,)

    @property
    def is_scalar(self):
        return self.lorentz_irrep == (0, 0)

    @property
    def is_boson(self):
        return sum(self.lorentz_irrep) % 2 == 0

    @property
    def is_left_fermion(self):
        return self.lorentz_irrep == (1, 0)

    @property
    def is_right_fermion(self):
        return self.lorentz_irrep == (0, 1)

    @property
    def is_fermion(self):
        return not self.is_boson

    @property
    def is_vector(self):
        return self.lorentz_irrep == (1, 1)

    @property
    def is_singlet(self):
        return self.dynkin == "00000" and self.y == 0

    @property
    def is_sm_singlet(self):
        return self.sm_irrep == (0, 0, 0) and self.y == 0

    @property
    def is_real_sm_irrep(self):
        is_real_colour = self.colour_irrep == tuple(reversed(self.colour_irrep))
        is_real_isospin = self.isospin_irrep[0] % 2 == 0
        is_real_y = self.y == 0
        return is_real_colour and is_real_isospin and is_real_y

    def __str__(self):
        return self.__repr__()

    @property
    def _dict(self):
        return {
            "label": self.label,
            "dynkin": self.dynkin,
            "charges": self.charges,
            "is_conj": self.is_conj,
            "comm": self.comm,
            "symmetry": self.symmetry,
            "nf": self.nf,
            "derivs": self.derivs,
            "stripped": self.stripped,
            # "history": self.history,
            # "latex": self.latex,
        }

    def __hash__(self):
        unhashable = "charges", "symmetry", "stripped"
        values = [v for k, v in self._dict.items() if not k in unhashable]

        # constructe hashable structures for each unhashable
        symmetry_tuple = tuple(map(tuple, self._dict["symmetry"]))
        charges_tuple = tuple(self.charges.items())

        data = tuple(values) + (symmetry_tuple,) + (charges_tuple,)
        return hash(data)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self._dict == other._dict

    def __mul__(self, other):
        grp = "SU2 x SU2 x SU3 x SU2"
        self_irrep = irrep(grp, " ".join(self.dynkin))
        other_irrep = irrep(grp, " ".join(other.dynkin))

        # basisgen returns a dict-like object: irrep -> multiplicity
        prod_dict = self_irrep * other_irrep

        # make sure charges are consistent
        assert self.charges.keys() == other.charges.keys()
        # add charges key-wise
        charges = {}
        for k, v in dict(self.charges).items():
            charges[k] = v + dict(other.charges)[k]

        # construct history
        new_hist = History(left=self, right=other)

        irreps = [
            Field(
                self.label_with_dagger + other.label_with_dagger,
                k.highest_weight.components,
                charges=charges,
                history=new_hist,
                multiplicity=v,
            )
            for k, v in prod_dict.items()
        ]
        return irreps

    @property
    def conj(self):
        dynkin = self.lorentz_irrep[::-1] + self.colour_irrep[::-1] + self.isospin_irrep
        return self.__class__(
            label=self.label,
            dynkin="".join(str(d) for d in dynkin),
            charges={k: -v for k, v in self.charges.items()},
            is_conj=(not self.is_conj),
            symmetry=self.symmetry,
            history=self.history,
            multiplicity=self.multiplicity,
            comm=self.comm,
            latex=self.latex,
            nf=self.nf,
            derivs=self.derivs,
            stripped=self.stripped,
        )

    def walked(self) -> Prod:
        if not self.history:
            return self

        return Prod(
            irrep=self,
            left=self.history.left.walked(),
            right=self.history.right.walked(),
        )

    def pprint(self) -> None:
        print(repr_tree(self.walked()))

    def fresh_indices(self) -> "IndexedField":
        # fresh_indices = " ".join(i for i in self.get_fresh_indices(store))
        fresh_indices = Index.fresh_indices(self.dynkin)
        label = (
            self.label_with_dagger if not isinstance(self, IndexedField) else self.label
        )

        g = ""
        if self.nf > 1:
            g += " "
            g += str(Index.fresh("g"))

        return IndexedField(
            label=label,
            indices=fresh_indices + g,
            charges=self.charges,
            is_conj=self.is_conj,
            symmetry=self.symmetry,
            comm=self.comm,
            latex=self.latex,
            nf=self.nf,
            derivs=self.derivs,
            stripped=self.stripped,
        )

    def get_latex(self):
        if self.is_conj:
            # add twidle for doublets so that their SU2 indices are always raised
            if self.isospin_irrep != (0,):
                if self.derivs == 0:
                    return rf"\tilde{{{self.latex}}}"
                # tilde will be on field get_latex method
                return self.latex

            # add dagger for singlets
            if self.derivs != 0:
                return self.latex
            return "{" + self.latex + r"^{\dagger}}"

        return self.latex

    @property
    def mass_dim(self):
        if self.is_fermion:
            return 1.5
        if self.is_scalar or self.is_vector:
            return 1

        raise Exception(
            f"Unknown mass dimension for lorentz irrep {self.lorentz_irrep}"
        )

    @classmethod
    def dynkin_difference(cls, dynkin_1: str, dynkin_2: str):
        """dynkin_2 - dynkin_1"""
        import numpy as np

        dyn_array_1 = np.array(list(map(int, dynkin_1)))
        dyn_array_2 = np.array(list(map(int, dynkin_2)))

        return list(dyn_array_2 - dyn_array_1)


class IndexedField(tensor.Tensor, Field):
    def __new__(cls, label: str, indices: str, symmetry=None, comm=0, **kwargs):
        if isinstance(indices, str):
            indices = indices.split()

        # classify index types and indices by name
        # e.g. 'i0' -> isospin, 'c0' -> colour, etc.
        tensor_indices = [Index(i) for i in indices]
        index_types = [i.tensor_index_type for i in tensor_indices]

        if symmetry is None:
            n_gauge_indices = (
                len(indices) if GENERATION not in index_types else len(indices) - 1
            )
            symmetry = [[1] * n_gauge_indices] if indices else []

        if isinstance(label, tensor.TensorHead):
            tensor_head = label
            label = str(label.args[0])
        else:
            sym = symmetry + ([[1]] if GENERATION in index_types else [])
            tensor_head = tensor.tensorhead(label, index_types, sym=sym, comm=comm)

        return super(IndexedField, cls).__new__(cls, tensor_head, tensor_indices)

    def __init__(
        self,
        label,
        indices,
        charges=None,
        is_conj=False,
        symmetry=None,
        comm=0,
        latex=None,
        nf=1,
        derivs=0,
        stripped=None,
        **kwargs,
    ):
        """Initialises IndexedField object.

        label: str (choice made to include the dagger)
        indices: space separated str or list of str
        charges: dict (must contain "y" as key)
        is_conj: bool

        Dynkin information from Field used to create correct indices for
        IndexedField.

        """
        # initialise charges again (in case initialised independently)
        # if charges is None:
        #     charges = {"y": 0, "3b": 0}

        # make sure charges contains hypercharge
        # assert "y" in charges

        if isinstance(label, tensor.TensorHead):
            label = str(label.args[0])

        Field.__init__(
            self,
            label=label,
            dynkin=get_dynkin(indices),
            charges=charges,
            is_conj=is_conj,
            symmetry=symmetry,
            comm=comm,
            latex=latex,
            nf=nf,
            derivs=derivs,
            stripped=stripped,
        )

        # self.comm = FERMI if self.is_fermion else BOSE
        self.index_labels = indices

    @property
    def field(self):
        return Field(
            label=self.label if not self.is_conj else self.label[:-1],
            dynkin=self.dynkin,
            charges=self.charges,
            is_conj=self.is_conj,
            symmetry=self.symmetry,
            comm=self.comm,
            latex=self.latex,
            nf=self.nf,
            derivs=self.derivs,
            stripped=self.stripped,
        )

    def substitute_indices(self, *indices):
        """Substitute indices according to ``indices`` with side effect."""
        new_indices = copy(self.indices)
        for i, j in indices:
            for pos, idx in enumerate(self.indices):
                if i == idx:
                    self.indices[pos] = j

    def substituted_indices(self, *indices):
        """Return a copy of the IndexedField with indices substituted according to ``indices``."""
        new_indices = copy(self.indices)
        for i, j in indices:
            for pos, idx in enumerate(self.indices):
                if i == idx:
                    new_indices[pos] = j

        return IndexedField(
            label=self.label,
            indices=new_indices,
            charges=self.charges,
            is_conj=self.is_conj,
            symmetry=self.symmetry,
            multiplicity=self.multiplicity,
            comm=self.comm,
            latex=self.latex,
            nf=self.nf,
            derivs=self.derivs,
            stripped=self.stripped,
        )

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
            latex=self.latex,
            comm=self.comm,
            nf=self.nf,
            derivs=self.derivs,
            stripped=self.stripped,
        )

    @property
    def indices_by_type(self):
        """Returns a dictionary mapping index type to tuple of indices.

        """
        # result = {k: [] for k in Index.get_index_types().values()}
        # for i in self.indices:
        #     result[i.index_type].append(i)
        # return {k: tuple(v) for k, v in result.items()}
        return Index.indices_by_type(indices=self.indices)

    @property
    def _dynkins(self):
        """Returns a dictionary of 2-tuples of integers mapping index type to (raised,
        lowered) indices.

        """
        return {k: Index.dynkin(v) for k, v in self.indices_by_type.items()}

    @property
    def _dict(self):
        return {
            "label": self.label,
            "indices": self.index_labels,
            "charges": self.charges,
            "is_conj": self.is_conj,
            "symmetry": self.symmetry,
            "comm": self.comm,
            "latex": self.latex,
            "nf": self.nf,
            "derivs": self.derivs,
            "stripped": self.stripped,
        }

    @property
    def pmatch_data(self):
        """A simple representation of the field for use with pattern matching.

        >>> pmatch_data(L("u0 i0"))
        ('L', ('u0', 'i0'), (('y', -1/2), ('3b', 0)))

        """
        return (
            self._dict["label"],
            tuple(self.indices_by_type.items()),
            tuple(sorted(self._dict["charges"].items(), key=lambda x: x[0])),
            self._dict["comm"],
            self.derivs,
            self.is_conj,
        )

    def __hash__(self):
        charges_tuple = tuple(self.charges.items())
        return hash(
            (
                self.label,
                tuple(self.indices),
                self.dynkin,
                charges_tuple,
                self.is_conj,
                self.comm,
            )
        )

    def __mul__(self, other):
        prod = normalise_operator_input(self, other)
        return Operator(*prod)

    def __add__(self, other):
        pass

    def __repr__(self):
        sympy_repr = super(IndexedField, self).__repr__()
        # if self.is_conj:
        #     split_sympy_repr = sympy_repr.split("(")
        #     new_repr = [split_sympy_repr[0] + "†("] + split_sympy_repr[1:]
        #     return "".join(new_repr)
        return sympy_repr

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False

        return self._dict == other._dict

    def __lt__(self, other):
        """An arbitrary way to order indexed fields."""
        comp_func = lambda x: f"{x.label}{x.dynkin}{x.indices}{x.charges}{x.is_conj}"
        return comp_func(self) < comp_func(other)

    def latex_and_pop(self, index_dict: dict, style: dict) -> str:
        """Returns the latex form of the indexed field while making state changes.

        Assigns unused indices of a certain type according to style. Pops new
        indices off style and saves the relation between Index objects and their
        representation in the index_dict (by side effect again).

        Keys in style dictate which indices are processed.

        """
        # Make sure latex attribute set
        if self.latex is None:
            raise ValueError(f"No latex string assigned to {self}")

        # deal with indices
        indices = []
        for index in self.indices:
            type_ = Index.get_index_labels()[index.index_type]
            # condition that style wants to be printed
            if type_ in style:
                index_string = style[type_].pop(0)
                index_dict[index] = index_string
                indices.append(index_string)

        if not indices:
            return self.get_latex()

        return rf"{self.get_latex()}^{{{' '.join(indices)}}}"

    def strip_derivs(self):
        if self.derivs == 0:
            return self

        other = {
            "is_conj": self.is_conj,
            "comm": self.comm,
            "nf": self.nf,
            "derivs": 0,
            "charges": {k: v for k, v in self.charges.items()},
            "stripped": self.stripped,
        }

        kwargs = {**self.stripped, **other}
        return Field(**kwargs)

    @property
    def gauge_indices(self):
        return [i for i in self.indices if i.index_type not in ("Dotted", "Undotted")]

    def strip_derivs_with_indices(self):
        """Remove the derivatives from a field and call on indices. New indices have a
        partial derivative symbol appended to the name.

        """
        if self.derivs == 0:
            return self

        field = self.strip_derivs()
        undotted, dotted, colour, isospin, _ = self.indices_by_type.values()
        indices = " ".join(str(i) for i in colour + isospin)
        lorentz = undotted + dotted
        if field.is_fermion and self.derivs % 2 == 0:
            lorentz = lorentz[0]
            indices = str(lorentz) + " " + indices
        elif field.is_fermion:
            lorentz = lorentz[0].conj
            # add partial symbol to index name to show that it came from derivative
            lorentz = str(lorentz) + "∂"
            indices = str(lorentz) + " " + indices

        return field(indices)


class Operator(tensor.TensMul):
    def __new__(cls, *args, **kwargs):
        return super(Operator, cls).__new__(cls, *args, **kwargs)

    def __init__(self, *args, **kwargs):
        self.tensors = [arg for arg in args if arg != 1]

    @property
    def _dict(self):
        return {"tensors": self.tensors}

    @property
    def contains_derivative(self):
        for t in self.tensors:
            if str(t).startswith("D"):
                return True
            if str(t).startswith("G"):
                return True
            if str(t).startswith("W"):
                return True
            if str(t).startswith("B"):
                return True

        return False

    @property
    def BL_classification(self):
        if self.contains_derivative:
            return -1

        sorted_fields = tuple(
            sorted(f.label_with_dagger.replace("†", ".conj") for f in self.fields)
        )
        if sorted_fields not in BL_LIST:
            return -1

        return BL_LIST[sorted_fields]

    @property
    def free_indices(self):
        return self.get_free_indices()

    # @property
    # def fields(self):
    #     return [f for f in self.args if not is_invariant_symbol(f)]

    @property
    def indexed_fields(self):
        return [f for f in self.tensors if isinstance(f, IndexedField)]

    @property
    def epsilons(self):
        return [f for f in self.tensors if not isinstance(f, IndexedField)]

    @property
    def fields(self):
        return [f.field for f in self.indexed_fields]

    def pickle_form(self):
        """Avoid trouble pickling objects caused by sympy inheritance and just pickle a
        state dictionary.

        """
        return {
            "fields": [f._dict for f in self.indexed_fields],
            "epsilons": [[i.label for i in e.indices] for e in self.epsilons],
        }

    @classmethod
    def from_pickle_form(cls, pickle_form: dict):
        """Avoid trouble unpickling objects caused by sympy inheritance and just read
        from a state dictionary.

        """
        fields = pickle_form["fields"]
        epsilons = [eps(" ".join(indices)) for indices in pickle_form["epsilons"]]
        indexed_fields = [IndexedField(**d) for d in fields]
        return cls(*indexed_fields, *epsilons)

    @property
    def structures(self):
        """Like epsilons but returns sympy tensor objects."""
        return [f for f in self.args if is_invariant_symbol(f)]

    def by_structures(self):
        order = Index.get_tensor_index_types().values()
        out = []
        for index_type in order:
            out.append([])
            for symbol in self.structures:
                if symbol.indices[0].tensor_index_type == index_type:
                    out[-1].append(symbol)

        return [l for l in out if l]

    def by_epsilon_structures(self):
        order = Index.get_tensor_index_types().values()
        out = []
        for index_type in order:
            out.append([])
            for symbol in self.epsilons:
                if symbol.indices[0].tensor_index_type == index_type:
                    out[-1].append(symbol)

        return [l for l in out if l]

    @property
    def duplicate_fields(self):
        """Returns a dictionary mapping the field label to the duplicates."""
        key_func = lambda f: str(f.args[0].args[0])
        fields = [t for t in self.args if not str(t).startswith("Eps")]
        sorted_fields = sorted(fields, key=key_func)
        out = {}
        for k, g in groupby(sorted_fields, key=key_func):
            g = list(g)
            if len(g) > 1:
                out[k] = flatten([f.indices for f in g])

        return out

    @property
    def safe_nocoeff(self):
        return self.nocoeff if not isinstance(self, Zero) else 0

    def simplify(self, fill=False):
        if fill:
            simple = self.fill_free_indices().sorted_components().canon_bp()
        else:
            simple = self.sorted_components().canon_bp()

        for label, tensor_index_type in Index.get_tensor_index_types().items():
            if label in {"u", "d", "i"}:
                if isinstance(simple, Zero):
                    return 0
                else:
                    simple = simple.contract_metric(tensor_index_type.metric)

        return simple.canon_bp()

    def safe_simplify(self):
        try:
            return self.simplify()
        except IndexError:
            # deal with no index case (sometimes throws an error)
            # collect fields with indices and only use canon_bp on those
            unindexed, indexed = [], []
            for f in self.tensors:
                if f.get_indices():
                    indexed.append(f)
                else:
                    unindexed.append(f)

            if not indexed:
                return self

            indexed_op = reduce(lambda x, y: x * y, indexed)
            return reduce(lambda x, y: x * y, unindexed, indexed_op.simplify())

    def is_equivalent_to(self, other: "Operator"):
        """Returns true if operators equal up to numerical constant or free index
        relabelling.

        """
        # Try naive thing for now
        return self.nocoeff == other.nocoeff

    def fill_free_indices(self):
        store = defaultdict(int)
        replacements = []
        for index_type, idxs in self.indices_by_type.items():
            for i in idxs:
                n = store[Index.get_index_labels()[index_type]]
                new_index = Index.cons_index(index_type, n, i.is_up)

                store[Index.get_index_labels()[index_type]] += 1
                replacements.append((i, new_index))

        tensors = [t.fun_eval(*replacements) for t in self.tensors]
        prod = normalise_operator_input(*tensors)
        return Operator(*prod)

    @property
    def indices_by_type(self):
        """Returns free indices in operator by type."""
        return Index.indices_by_type(indices=self.free_indices)

    @property
    # Will only return something meaninful if no extra epsilons or deltas
    # around. TODO Extend to ignore uncontracted epsilons and deltas (from
    # is_contracted_epsilon in completions.completions)
    def dynkin(self):
        return get_dynkin(" ".join(str(i) for i in self.free_indices))

    @property
    def dynkin_ints(self):
        return [int(i) for i in self.dynkin]

    # Always multiply invariant symbols on the right
    def __mul__(self, other):
        prod = normalise_operator_input(*self.tensors, other)
        return Operator(*prod)

    def latex(self, ignore=None):
        # ordering of fields within the oeprator and style of indices
        field_ordering = ["L", "e", "Q", "u", "d", "H", "B", "W", "G"]
        for field in self.indexed_fields:
            f, fc = field, field.conj
            if hasattr(f, "label") and not f.label in field_ordering:
                field_ordering.append(f.label.replace("~", ""))
                field_ordering.append(fc.label.replace("~", ""))

        # indices i, j, ... q used for isospin
        isospin_indices = list(ascii_lowercase[8:])
        isospin_indices.remove("o")  # remove o because it's unsightly

        style = {
            "i": isospin_indices,
            "c": list(ascii_lowercase[:8]),
            "u": copy(TEX_GREEK_LOWERCASE),
            "d": copy(DOTTED_TEX_GREEK_LOWERCASE),
            # "g": list(ascii_lowercase[19:]),
        }

        # remove styles you don't care about
        if ignore is not None:
            for char in ignore:
                style.pop(char)

        # extract first non D in field name e.g. H for DDH
        first_non_deriv = lambda f: str(f).split("D")[-1][0]
        order_func = lambda f: field_ordering.index(first_non_deriv(f))
        sorted_fields = sorted(self.indexed_fields, key=order_func)

        # maps Index objects to string of latex index
        index_dict = {}

        latex_strings = []
        for indexed_field in sorted_fields:
            field_latex = indexed_field.latex_and_pop(index_dict, style)
            latex_strings.append(field_latex)

        latex_epsilons = []
        for epsilon in self.epsilons:
            symb = r"\epsilon"
            # tread deltas differently
            if str(epsilon)[0] == "K":
                symb = r"\delta"

            # If two epsilons contracted, they will have an additional index
            # that isn't in index_dict, need to add it in manually. For now,
            # disregard raised and lowered indices for SU(3)
            eps_indices = []
            for idx in epsilon.indices:
                if -idx not in index_dict.keys():
                    type_ = Index.get_index_labels()[idx.index_type]
                    index_string = style[type_].pop(0)
                    index_dict[-idx] = index_string
                    index_dict[idx] = index_string

                eps_indices.append(index_dict[-idx])

            eps_indices.sort()
            # eps_indices = sorted([index_dict[-idx] for idx in epsilon.indices])
            eps_latex = rf"{symb}_{{{' '.join(eps_indices)}}}"
            latex_epsilons.append(eps_latex)

        # add on epsilons
        if latex_epsilons:
            latex_strings += [r" \cdot "] + sorted(latex_epsilons)

        return " ".join(latex_strings)

    def _repr_html_(self):
        return f"${self.latex()}$"

    def __deepcopy__(self, memo):
        return self.__class__(tensors=deepcopy(self.tensors, memo))

    @property
    def mass_dim(self):
        return sum(f.mass_dim for f in self.fields)


def assert_consistent_indices(
    indices: List[str],
    index_types: List[tensor.TensorIndexType],
    colour: Tuple[int, int],
):
    assert len(indices) == len(index_types)
    up_count, down_count = 0, 0
    for i, t in zip(indices, index_types):
        # make sure index names match
        assert Index(i).tensor_index_type == t

        # keep track of raised and lowered indices for colour
        if t == COLOUR:
            if i[0] == "-":
                down_count += 1
            else:
                up_count += 1
        else:
            # make sure no lowered indices for the SU2s
            assert i[0] != "-"

    assert (up_count, down_count) == colour


def get_dynkin(indices: Union[str, Iterable[Index]]):
    """get_dynkin("u0 c0 -c1 i0") => '10111'"""
    dynkin = {
        "Undotted": {True: 0},
        "Dotted": {True: 0},
        "Colour": {True: 0, False: 0},
        "Isospin": {True: 0},
    }

    if isinstance(indices, str):
        indices = indices.split()

    for i in indices:
        idx = Index(i)

        if idx.index_type == "Generation":
            continue

        dynkin[idx.index_type][idx.is_up] += 1

    flat_dynkin = flatten(map(lambda x: list(x.values()), dynkin.values()))
    return "".join(str(x) for x in flat_dynkin)


# @lru_cache(maxsize=None)
def decompose_product(*fields) -> List[Field]:
    """Decompose product of Fields.

    Example:
        >>> decompose_product(A, B, A.conj)
        >>> [ABA†(23003)(1/6), ABA†(23001)(1/6), ABA†(21003)(1/6), ...]

    """

    if len(fields) == 1:
        return fields[0]

    if len(fields) == 2:
        fst, snd = fields
        return fst * snd

    fst, snd, *rst = fields
    result = map(lambda x: decompose_product(x, *rst), fst * snd)
    return flatten(result)


def eps(indices: str):
    indices = [Index(i) for i in indices.split()]
    # make sure same index type
    t = indices[0].index_type
    tensor_index_type = indices[0].tensor_index_type

    # ensure valid index type
    assert t in Index.get_index_types().values()

    # return metric for the su2s but epsilon for su3
    if not t == Index.get_index_types()["c"]:
        # check consistent SU(2) indices
        assert len(indices) == 2
        for i in indices:
            assert not i.is_up

        return tensor_index_type.metric(*indices)

    # check consistent SU(3) indices
    assert len(indices) == 3
    position = indices[0].is_up
    for i in indices:
        assert position == i.is_up
    return tensor_index_type.epsilon(*indices)


def delta(indices: str):
    indices = [Index(i) for i in indices.split()]

    # check consistent SU(3) delta
    assert len(indices) == 2
    i, j = indices
    assert i.index_type == Index.get_index_types()["c"]
    assert j.index_type == Index.get_index_types()["c"]
    assert i.is_up
    assert not j.is_up

    return Index.get_tensor_index_types()["c"].delta(*indices)


def is_invariant_symbol(tensor):
    str_tensor = str(tensor)
    return (
        str_tensor.startswith("Eps")
        or str_tensor.startswith("KD")
        or str_tensor.startswith("metric")
    )


def normalise_operator_input(*args):
    """Input to Operator may contain products of epsilons which are sympy TensMul
    objects. This causes problems since I don't want to define products of
    epsilons or deltas to be operators.

    This function normalises the input to the Operator constructor to avoid
    problems.

    """
    new_product = []
    for fac in args:
        if isinstance(fac, tensor.TensMul):
            new_product += list(fac.args)
        else:
            new_product.append(fac)

    return sorted(new_product, key=type)


def D(field, dynkin):
    """A derivative.

    Returns a new field with additional dotted and undotted indices.

    Example:
       >>> D(L, "01")
       DL(01001)(-1/2)

       >>> D(L, "21")
       DL(21001)(-1/2)

    """
    undotted_delta = int(dynkin[0]) - field.dynkin_ints[0]
    dotted_delta = int(dynkin[1]) - field.dynkin_ints[1]

    # derivative can only change one dotted and one undotted index
    assert abs(undotted_delta) == 1
    assert abs(dotted_delta) == 1

    # other info to construct field instance
    deriv_symbol = "D"
    symbol = deriv_symbol + field.label
    new_field_dynkin = dynkin + field.dynkin[2:]
    rest = {
        "charges": field.charges,
        "comm": field.comm,
        "is_conj": field.is_conj,
        "nf": field.nf,
        "stripped": field.stripped,
    }

    new_field = Field(symbol, dynkin=new_field_dynkin, **rest)
    new_field.latex = f"(D{strip_parens(field.get_latex())})"
    new_field.derivs = field.derivs + 1
    # only add this information for the first derivative
    if new_field.stripped is None:
        new_field.stripped = {
            "label": field.label,
            "dynkin": field.dynkin,
            "symmetry": field.symmetry,
            "charges": field.charges,
            "latex": field.latex,
        }
    return new_field
