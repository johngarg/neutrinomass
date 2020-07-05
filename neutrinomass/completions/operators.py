#!/usr/bin/env python3

import pickle
import os

from functools import reduce

from neutrinomass.tensormethod.core import eps, Operator
from neutrinomass.tensormethod.sm import L, Q, db, ub, eb, H

from neutrinomass.completions.core import EffectiveOperator


def prod(lst):
    return reduce(lambda x, y: x * y, lst)


# read in pickled data from tensormethod script lnvlatex
with open(os.path.join(os.path.dirname(__file__), "operators.p"), "rb") as f:
    pickle_form_ops = pickle.load(f)

with open(os.path.join(os.path.dirname(__file__), "deriv_operators.p"), "rb") as f:
    pickle_form_deriv_ops = pickle.load(f)

# define operators
EFF_OPERATORS = {}
for k, op in pickle_form_ops.items():
    EFF_OPERATORS[k] = EffectiveOperator(k, Operator.from_pickle_form(op))

DERIV_EFF_OPERATORS = {}
for k, op in pickle_form_deriv_ops.items():
    if k in ("D1", "D11"):  # non-explosive
        continue
    if k.startswith("D19"):  # four-deriv
        continue

    DERIV_EFF_OPERATORS[k] = EffectiveOperator(k, Operator.from_pickle_form(op))

# add additional operators not in BL/dGJ list

# aux. Operator objects
o1 = EFF_OPERATORS["1"].operator.tensors
o2 = EFF_OPERATORS["2"].operator.tensors
o3a = EFF_OPERATORS["3a"].operator.tensors
o3b = EFF_OPERATORS["3b"].operator.tensors
o4a = EFF_OPERATORS["4a"].operator.tensors
o4b = EFF_OPERATORS["4b"].operator.tensors
o5a = EFF_OPERATORS["5a"].operator.tensors
o5b = EFF_OPERATORS["5b"].operator.tensors
o6a = EFF_OPERATORS["6a"].operator.tensors
o6b = EFF_OPERATORS["6b"].operator.tensors
o7 = EFF_OPERATORS["7"].operator.tensors
o8 = EFF_OPERATORS["8"].operator.tensors
o61a = EFF_OPERATORS["61a"].operator.tensors
o71 = EFF_OPERATORS["71"].operator.tensors
o76 = [
    eb.conj("d0"),
    eb.conj("d1"),
    ub.conj("d2 c0"),
    ub.conj("d3 c1"),
    db("u0 -c2"),
    db("u1 -c3"),
]
o82 = [
    L("u0 i0"),
    L.conj("d0 i1"),
    eb.conj("d1"),
    eb.conj("d2"),
    ub.conj("d3 c0"),
    db("u1 -c1"),
    H("i2"),
    H("i3"),
    eps("-i0 -i2"),
    eps("-i1 -i3"),
]
oyec = [L.conj("d3 i6"), eb.conj("d4"), H("i7"), eps("-i6 -i7")]
oydc = [Q.conj("d3 -c0 i6"), db.conj("d4 c1"), H("i7"), eps("-i6 -i7")]
prime = [H("i4"), H.conj("i5"), eps("-i4 -i5")]
prime_prime = [
    H("i4"),
    H.conj("i5"),
    H("i6"),
    H.conj("i7"),
    eps("-i6 -i7"),
    eps("-i4 -i5"),
]
prime_prime_prime = [
    H("i4"),
    H.conj("i5"),
    H("i6"),
    H.conj("i7"),
    H("i8"),
    H.conj("i9"),
    eps("-i8 -i9"),
    eps("-i6 -i7"),
    eps("-i4 -i5"),
]

# new operators from table in paper
EFF_OPERATORS["77"] = EffectiveOperator("77", prod(o1 + oyec))
EFF_OPERATORS["78"] = EffectiveOperator("78", prod(o1 + oydc))
EFF_OPERATORS["1p"] = EffectiveOperator("1p", prod(o1 + prime))
EFF_OPERATORS["8p"] = EffectiveOperator("8p", prod(o8 + prime))
EFF_OPERATORS["1pp"] = EffectiveOperator("1pp", prod(o1 + prime_prime))
# EFF_OPERATORS["1ppp"] = EffectiveOperator("1ppp", prod(o1 + prime_prime_prime))
EFF_OPERATORS["7p"] = EffectiveOperator("7p", prod(o7 + prime))
EFF_OPERATORS["8pp"] = EffectiveOperator("8pp", prod(o8 + prime_prime))
EFF_OPERATORS["71p"] = EffectiveOperator("71p", prod(o71 + prime))
EFF_OPERATORS["76p"] = EffectiveOperator("76p", prod(o76 + prime))
EFF_OPERATORS["77p"] = EffectiveOperator("77p", prod(o1 + oyec + prime))
EFF_OPERATORS["78p"] = EffectiveOperator("78p", prod(o1 + oydc + prime))
EFF_OPERATORS["79a"] = EffectiveOperator("79a", prod(o61a + prime))
EFF_OPERATORS["79b"] = EffectiveOperator("79b", prod(o2 + prime_prime))
EFF_OPERATORS["80a"] = EffectiveOperator("80a", prod(o5a + prime))
EFF_OPERATORS["80b"] = EffectiveOperator("80b", prod(o5b + prime))
EFF_OPERATORS["80c"] = EffectiveOperator("80c", prod(o3a + prime_prime))
EFF_OPERATORS["80d"] = EffectiveOperator("80d", prod(o3b + prime_prime))
EFF_OPERATORS["81a"] = EffectiveOperator("81a", prod(o6a + prime))
EFF_OPERATORS["81b"] = EffectiveOperator("81b", prod(o6b + prime))
EFF_OPERATORS["81c"] = EffectiveOperator("81c", prod(o4a + prime_prime))
EFF_OPERATORS["81d"] = EffectiveOperator("81d", prod(o4b + prime_prime))
EFF_OPERATORS["82"] = EffectiveOperator("82", prod(o82))
