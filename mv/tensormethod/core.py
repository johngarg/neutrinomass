#!/usr/bin/env python

import sympy.tensor.tensor as tensor

ISOSPIN = tensor.TensorIndexType("Isospin", metric=True, dummy_fmt="I", dim=2)
COLOUR = tensor.TensorIndexType("Colour", metric=None, dummy_fmt="C", dim=3)
GENERATION = tensor.TensorIndexType("Generation", metric=None, dummy_fmt="G", dim=3)
UNDOTTED = tensor.TensorIndexType("Undotted", metric=True, dummy_fmt="U", dim=2)
DOTTED = tensor.TensorIndexType("Dotted", metric=True, dummy_fmt="D", dim=2)


class Index(tensor.TensorIndex):
    def __new__(cls, label: str, *args, **kwargs):
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
        self.label = label  # the label that was passed in

    def __neg__(self):
        is_up = self.is_up
        new_label = "-" + self.label if is_up else self.label[1:]
        return Index(label=new_label)

    @property
    def conj(self):
        if self.index_type in Index.get_lorentz_index_types().values():
            if self.index_type == "Undotted":
                return Index(label=self.label.replace("u", "d"))
            return Index(label=self.label.replace("d", "u"))
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


class Field(tensor.Tensor):
    def __new__(cls, label: str, indices: str, symmetry=None, **kwargs):
        if isinstance(indices, str):
            indices = indices.split()

        if symmetry is None:
            symmetry = [[1]] * len(indices)

        # classify index types and indices by name
        # e.g. 'i0' -> isospin, 'c0' -> colour, etc.
        tensor_indices = [Index(i) for i in indices]
        index_types = [i.tensor_index_type for i in tensor_indices]
        tensor_head = tensor.tensorhead(label, index_types, sym=symmetry)

        return super(Field, cls).__new__(cls, tensor_head, tensor_indices)

    def __init__(self, label, indices, charges=None, is_conj=False, **kwargs):
        """Initialises Field object.

        label: str
        indices: space separated str or list of str
        charges: dict (must contain "y" as key)
        is_conj: bool

        """
        # Initialise field with only hypercharge (= to 0)
        if charges is None:
            charges = {"y": 0}
        if "y" not in charges.keys():
            charges = {**charges, "y": 0}

        self.label = label
        self.index_names = indices.split() if isinstance(indices, str) else indices
        self.charges = charges
        self.y = charges["y"]  # hypercharge
        self.is_conj = is_conj

    @property
    def conj(self):
        """Returns a copy of self but conjugated"""
        return Field(
            label=self.label,
            indices=[i.conj.label for i in self.indices],
            charges={k: -v for k, v in self.charges.items()},
            is_conj=(not self.is_conj),
        )

    @property
    def indices_by_type(self):
        result = {k: [] for k in Index.get_index_types().values()}
        for i in self.indices:
            result[i.index_type].append(i)
        return {k: tuple(v) for k, v in result.items()}

    @property
    def dynkins(self):
        return {k: Index.dynkin(v) for k, v in self.indices_by_type.items()}

    @property
    def sm_dynkins(self):
        dynkins = self.dynkins
        return dynkins["Colour"], dynkins["Isospin"]

    @property
    def lorentz_dynkins(self):
        dynkins = self.dynkins
        return dynkins["Undotted"], dynkins["Dotted"]

    @property
    def quantum_numbers(self):
        return self.sm_dynkins + (self.y,)

    @property
    def is_scalar(self):
        return self.lorentz_dynkins == (0, 0)

    @property
    def is_left_fermion(self):
        return self.lorentz_dynkins == (1, 0)

    @property
    def is_right_fermion(self):
        return self.lorentz_dynkins == (0, 1)

    @property
    def is_fermion(self):
        return self.is_left_fermion or self.is_right_fermion

    @property
    def _dict(self):
        return (self.label, self.indices, self.charges, self.is_conj)

    def __mul__(self, other):
        pass

    def __add__(self, other):
        pass

    def __repr__(self):
        sympy_repr = super(Field, self).__repr__()
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
    pass
