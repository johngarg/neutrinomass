#!/usr/bin/env python


from collections import namedtuple
from functools import lru_cache
from typing import Iterable, List, NamedTuple, Tuple, Union

import sympy.tensor.tensor as tensor
from basisgen import irrep
from sympy import Rational, flatten

from utils import repr_tree


class History(NamedTuple):
    left: "Field"
    right: "Field"


class Prod(NamedTuple):
    irrep: "Field"
    left: Union["Prod", "Field"]
    right: Union["Prod", "Field"]


ISOSPIN = tensor.TensorIndexType("Isospin", metric=True, dummy_fmt="I", dim=2)
COLOUR = tensor.TensorIndexType("Colour", metric=None, dummy_fmt="C", dim=3)
GENERATION = tensor.TensorIndexType("Generation", metric=None, dummy_fmt="G", dim=3)
UNDOTTED = tensor.TensorIndexType("Undotted", metric=True, dummy_fmt="U", dim=2)
DOTTED = tensor.TensorIndexType("Dotted", metric=True, dummy_fmt="D", dim=2)


class Index(tensor.TensorIndex):
    def __new__(cls, label: str, *args, **kwargs):
        if isinstance(label, Index):
            return label

        # This will be called in __mul__ method of IndexedTensor
        # i.e. Will be a contracted index
        if isinstance(label, tensor.TensorIndex):
            label = str(label).lower()

        # deal with negative sign for down index
        is_up = True
        if label[0] == "-":
            is_up = False
            label = label[1:]

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
        if self.index_type == "Isospin":
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
    def get_tensor_index_types(cls):
        return {"u": UNDOTTED, "d": DOTTED, "c": COLOUR, "i": ISOSPIN, "g": GENERATION}

    @classmethod
    def get_index_types(cls):
        return {
            "u": "Undotted",
            "d": "Dotted",
            "c": "Colour",
            "i": "Isospin",
            "g": "Generation",
        }

    @classmethod
    def get_sm_index_types(cls):
        return {"c": "Colour", "i": "Isospin"}

    @classmethod
    def get_lorentz_index_types(cls):
        return {"u": "Undotted", "d": "Dotted"}

    @classmethod
    def classify_index(cls, idx: str) -> tensor.TensorIndexType:
        return Index.get_tensor_index_types()[idx[0]]


class Field:
    def __init__(
        self,
        label: str,
        dynkin: str,
        charges=None,
        is_conj=False,
        comm=0,
        symmetry=None,
        history=None,
        multiplicity=1,
        **kwargs,
    ):
        """Field("A", dynkin="1 0 0 1 1", charges={"y": 1})

        Tensor symmetry options:
            ``[[1]]``         vector
            ``[[1]*n]``       symmetric tensor of rank ``n``
            ``[[n]]``         antisymmetric tensor of rank ``n``
            ``[[2, 2]]``      monoterm slot symmetry of the Riemann tensor
            ``[[1],[1]]``     vector*vector
            ``[[2],[1],[1]]`` (antisymmetric tensor)*vector*vector

        """

        # initialise charges
        if charges is None:
            charges = {"y": 0}

        # make sure charges contains hypercharge
        assert "y" in charges.keys()

        # for SU(2) and SU(3), indices will always be symmetric
        if symmetry is None:
            symmetry = []
            for i in map(int, dynkin):
                if i:
                    symmetry += [[1] * i]

        if not isinstance(dynkin, str):
            dynkin = "".join(str(i) for i in dynkin)

        self.symmetry = symmetry
        self.history = history
        self.charges = charges
        self.multiplicity = multiplicity
        self.label = label
        self.dynkin = dynkin.strip()
        self.comm = comm
        self.is_conj = is_conj

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
        assert_consistent_indices(indices.split(), index_types, (su3_up, su3_down))
        return IndexedField(
            label=self.label,
            indices=indices,
            charges=self.charges,
            is_conj=self.is_conj,
            symmetry=self.symmetry,
            multiplicity=self.multiplicity,
            comm=self.comm,
        )

    @property
    def label_with_dagger(self):
        maybe_conj = "†" if self.is_conj else ""
        return self.label + maybe_conj

    def __repr__(self):
        return self.label_with_dagger + f"({self.dynkin})({self.y})"

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
    def is_left_fermion(self):
        return self.lorentz_irrep == (1, 0)

    @property
    def is_right_fermion(self):
        return self.lorentz_irrep == (0, 1)

    @property
    def is_fermion(self):
        return self.is_left_fermion or self.is_right_fermion

    @property
    def is_vector(self):
        return self.lorentz_irrep == (1, 1)

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
            # self.multiplicity,
            # self.history,
        }

    def __hash__(self):
        unhashable = "charges", "symmetry"
        values = [v for k, v in self._dict.items() if not k in unhashable]

        # constructe hashable structures for each unhashable
        symmetry_tuple = tuple(map(tuple, self._dict["symmetry"]))
        charges_tuple = tuple(self.charges.items())

        return hash(tuple(values) + (symmetry_tuple,) + (charges_tuple,))

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
        )

    @property
    def walked(self) -> Prod:
        if not self.history:
            return self

        return Prod(
            irrep=self, left=self.history.left.walked, right=self.history.right.walked
        )

    @property
    def pprint(self) -> None:
        print(repr_tree(self.walked))

    def get_fresh_indices(self) -> List[str]:
        # names of index types
        # don't include generation
        u, d, c, i = list(Index.get_index_types().keys())[:-1]
        # number of each index type
        nu, nd, nc_up, nc_down, ni = self.dynkin_ints

        indices = {}
        # lorentz
        indices[u] = [u + str(x) for x in store[:nu]]
        indices[d] = [d + str(x) for x in store[:nd]]
        # colour
        indices[c] = [c + str(x) for x in store[:nc_up]]
        indices[c] += ["-" + c + str(x) for x in store[nc_up : nc_up + nc_down]]
        # isospin
        indices[i] = [i + str(x) for x in store[:ni]]

        # make the choice to return flattened values
        # full dict available if you need it (just return indices below)
        return flatten(indices.values())

    def fresh_indexed_field(self, store: Iterable[int] = range(10)) -> "IndexedField":
        fresh_indices = " ".join(i for i in self.get_fresh_indices(store))
        return IndexedField(
            label=self.label,
            indices=fresh_indices,
            charges=self.charges,
            is_conj=self.is_conj,
            symmetry=self.symmetry,
            comm=self.comm,
        )


class IndexedField(tensor.Tensor, Field):
    def __new__(cls, label: str, indices: str, symmetry=None, **kwargs):
        if isinstance(indices, str):
            indices = indices.split()

        if symmetry is None:
            symmetry = [[1] * len(indices)]

        # classify index types and indices by name
        # e.g. 'i0' -> isospin, 'c0' -> colour, etc.
        tensor_indices = [Index(i) for i in indices]
        index_types = [i.tensor_index_type for i in tensor_indices]
        if isinstance(label, tensor.TensorHead):
            tensor_head = label
        else:
            tensor_head = tensor.tensorhead(label, index_types, sym=symmetry)

        return super(IndexedField, cls).__new__(cls, tensor_head, tensor_indices)

    def __init__(
        self,
        label,
        indices,
        charges=None,
        is_conj=False,
        symmetry=None,
        comm=0,
        **kwargs,
    ):
        """Initialises IndexedField object.

        label: str
        indices: space separated str or list of str
        charges: dict (must contain "y" as key)
        is_conj: bool

        Dynkin information from Field used to create correct indices for
        IndexedField.

        """
        # initialise charges again (in case initialised independently)
        if charges is None:
            charges = {"y": 0}

        # make sure charges contains hypercharge
        assert "y" in charges

        Field.__init__(
            self,
            label=label,
            dynkin=get_dynkin(indices),
            charges=charges,
            is_conj=is_conj,
            symmetry=symmetry,
            comm=comm,
        )

        self.index_labels = indices

    @property
    def field(self):
        return Field(
            label=self.label,
            dynkin=self.dynkin,
            charges=self.charges,
            is_conj=self.is_conj,
            symmetry=self.symmetry,
            comm=self.comm,
        )

    @property
    def conj(self):
        """Returns a copy of self but conjugated"""
        return self.__class__(
            label=self.label,
            indices=" ".join(i.conj.label for i in self.indices),
            charges={k: -v for k, v in self.charges.items()},
            is_conj=(not self.is_conj),
            symmetry=self.symmetry,
            comm=self.comm,
        )

    @property
    def indices_by_type(self):
        """Returns a dictionary mapping index type to tuple of indices.

        """
        result = {k: [] for k in Index.get_index_types().values()}
        for i in self.indices:
            result[i.index_type].append(i)
        return {k: tuple(v) for k, v in result.items()}

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
            "index_labels": self.index_labels,
            "charges": self.charges,
            "is_conj": self.is_conj,
            "symmetry": self.symmetry,
            "comm": self.comm,
        }

    def __hash__(self):
        return super(self.__class__, self).__hash__()

    def __mul__(self, other):
        return Operator(self, other)

    def __add__(self, other):
        pass

    def __repr__(self):
        sympy_repr = super(self.__class__, self).__repr__()
        if self.is_conj:
            split_sympy_repr = sympy_repr.split("(")
            new_repr = [split_sympy_repr[0] + "†("] + split_sympy_repr[1:]
            return "".join(new_repr)
        return sympy_repr

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self._dict == other._dict


class Operator(tensor.TensMul):
    def __new__(cls, *args, **kwargs):
        return super(Operator, cls).__new__(cls, *args, **kwargs)

    def __init__(self, *args, **kwargs):
        self.tensors = args

    @property
    def free_indices(self):
        return self.get_free_indices()

    @property
    def dynkin(self):
        return get_dynkin("".join(str(i) for i in self.free_indices))


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
        dynkin[idx.index_type][idx.is_up] += 1

    flat_dynkin = flatten(map(lambda x: list(x.values()), dynkin.values()))
    return "".join(str(x) for x in flat_dynkin)


@lru_cache(maxsize=None)
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

    if not t == Index.get_index_types()["c"]:
        # check consistent SU(2) indices
        assert len(indices) == 2
        for i in indices:
            assert not i.is_up
    else:
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
