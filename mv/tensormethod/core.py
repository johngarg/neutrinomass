#!/usr/bin/env python

import sympy.tensor.tensor as tensor

ISOSPIN = tensor.TensorIndexType("Isospin", metric=True, dummy_fmt="I", dim=2)
COLOUR = tensor.TensorIndexType("Colour", metric=None, dummy_fmt="C", dim=3)
GENERATION = tensor.TensorIndexType("Generation", metric=None, dummy_fmt="G", dim=3)
UNDOTTED = tensor.TensorIndexType("Undotted", metric=True, dummy_fmt="U", dim=2)
DOTTED = tensor.TensorIndexType("Dotted", metric=True, dummy_fmt="D", dim=2)

INDEX_DICT = {"u": UNDOTTED, "d": DOTTED, "c": COLOUR, "i": ISOSPIN}


def classify_index(idx: str) -> tensor.TensorIndexType:
    if isinstance(idx, str):
        return INDEX_DICT[idx[0]]
    return idx.tensor_index_type


class Field(tensor.Tensor):
    def __new__(cls, label: str, indices: str, symmetry=None, **kwargs):
        if isinstance(indices, str):
            indices = indices.split()
        elif isinstance(indices, tuple):
            # called by sympy's Tensor methods
            label = label.name
            indices = [str(i) for i in indices]

        if symmetry is None:
            symmetry = [[1]] * len(indices)

        # classify index types and indices by name
        # e.g. 'i0' -> isospin, 'c0' -> colour, etc.
        index_types = [classify_index(i) for i in indices]
        tensor_indices = [
            # need to pass first argument as string because when called by
            # sympy, TensorIndex passed in
            tensor.tensor_indices(str(i), classify_index(i))
            for i in indices
        ]
        tensor_head = tensor.tensorhead(label, index_types, sym=symmetry)

        return super(Field, cls).__new__(cls, tensor_head, tensor_indices)

    def __init__(self, label, indices, y=0, **kwargs):
        self.label = label
        self.index_names = indices
        self.y = y  # hypercharge

    def __mul__(self, other):
        pass

    def __add__(self, other):
        pass


class Operator(tensor.TensMul):
    pass
