#!/usr/bin/env python3

from neutrinomass.tensormethod import D, L, Q, H, eb, ub, db, eps, delta
from neutrinomass.tensormethod.core import IndexedField, Field
from neutrinomass.completions.topologies import Leaf
from neutrinomass.completions.core import (
    EffectiveOperator,
    Completion,
    cons_completion_field,
)
from neutrinomass.utils.functions import conjugate_field, conjugate_term
from neutrinomass.database.utils import subsets
from functools import reduce
import operator

from networkx import Graph
from collections import defaultdict
import time
from typing import Dict, List
import re
import math
import numpy as np
from sympy import prime
from sympy.ntheory import factorint
from itertools import groupby
import os
from glob import glob
import pandas as pd
import pickle

ExoticField = cons_completion_field


def match(pattern: str, data: str):
    return re.match(pattern, data)


def matches_in(pattern, coll):
    for x in coll:
        if match(pattern, x):
            return True
    return False


def Partition(*args):
    return tuple(args)


class LazyCompletion:
    def __init__(self, head: dict, tail: str):
        """head is a dict with essential information about completion. For complete
        completion, call `force` method.

        """
        self.head = head
        self.tail = tail

        assert "quantum_numbers" in self.head
        assert "operator_name" in self.head
        # assert "operator_dimension" in self.head

    def force(self):
        return eval(self.tail)

    @property
    def quantum_numbers(self):
        return set(self.head["quantum_numbers"])

    @property
    def operator_name(self):
        return self.head["operator_name"]

    def contains_field(self, pattern):
        return matches_in(pattern, self.quantum_numbers)

    def contains_interaction(self, patterns):
        """This probably needs to be rewritten recursively"""
        interactions = self.head["terms"]
        for interaction in interactions:
            is_match = True
            for pattern in patterns:
                is_match = is_match and bool(matches_in(pattern, interaction))

            if is_match:
                return True

        return False

    # @property
    # def operator_dimension(self):
    #     return self.head["operator_dimension"]


def read_completions(filename: str):
    """Do this as a context manager?"""
    completions = defaultdict(list)
    with open(filename, "r") as f:
        line = f.readline()
        while line:
            comp = eval(line)
            completions[comp.operator_name].append(comp)
            line = f.readline()

    return completions


class ModelDataFrame(pd.DataFrame):
    _metadata = [
        "exotics",
        "terms",
        "exotic2int",
        "int2exotic",
        "term2int",
        "int2term",
    ]

    @property
    def _constructor(self):
        return ModelDataFrame

    @classmethod
    def new(cls, data, exotics, terms):
        df = cls(data)

        df.exotic2int = {k: v for k, v in exotics.items()}
        df.int2exotic = {v: k for k, v in exotics.items()}
        df.exotics = {**df.exotic2int, **df.int2exotic}

        df.term2int = {k: v for k, v in terms.items()}
        df.int2term = {v: k for k, v in terms.items()}
        df.terms = {**df.term2int, **df.int2term}

        return df

    def completion(self, index: int) -> Completion:
        return eval(self["completion"][index])

    def related_models(self, index: int) -> "ModelDataFrame":
        demo_num = self["democratic_num"][index]
        factor_groups = [l for l in subsets(factorint(demo_num)) if len(l) > 1]

        indices = []
        for group in factor_groups:
            prod = reduce(operator.mul, group)
            indices += list(self[self["democratic_num"] == prod].index)

        indices.remove(index)
        return self.loc[indices]


class ModelDatabase:
    def __init__(
        self,
        path: str,
        philosophy: str = "democratic",
        criterion: str = "mass",
        data=None,
    ):
        self.philosophy = philosophy
        self.criterion = criterion
        self.is_forced = False

        assert self.philosophy in {"democratic", "stringent"}
        assert self.criterion in {"mass", "dimension"}

        if data is None:
            self.path = path

            # print("Initialising database...")
            filenames = glob(os.path.join(self.path, "*.dat"))
            mvdb_data = [dict(read_completions(f)) for f in filenames]
            mv_dict = {k: v for d in mvdb_data for k, v in d.items()}
            self.data = mv_dict

        else:
            self.path = None
            self.is_ordered = True
            self.data = data

        # initialise prime dictionary
        exotics = {}
        counter = 1
        for k, v in self.data.items():
            for model in v:
                for exotic in model.quantum_numbers:
                    if exotic not in exotics:
                        exotics[conjugate_field(exotic)] = prime(counter)
                        exotics[exotic] = prime(counter)
                        counter += 1

        # dict mapping exotic field (str representation) to unique prime number
        self.exotic_prime_dict = exotics
        self.inv_exotic_prime_dict = {v: k for k, v in self.exotic_prime_dict.items()}

        term_dict = {}
        counter = 1
        for k, v in self.data.items():
            for model in v:
                n_terms = len(model.head["terms"])
                for i in range(n_terms):
                    # sort all of the terms by side effect
                    term = tuple(sorted(model.head["terms"][i]))
                    model.head["terms"][i] = term
                    if term not in term_dict:
                        # add conj term first so that when filtering you keep
                        # the unconjugated term
                        term_dict[conjugate_term(term)] = prime(counter)
                        term_dict[term] = prime(counter)
                        counter += 1

        # dict mapping sorted tuple of strings representing interaction (and
        # conjugate) to unique prime number
        self.term_prime_dict = term_dict
        self.inv_term_prime_dict = {v: k for k, v in self.term_prime_dict.items()}

        # dict mapping operator label to neutrino-mass scale estimate
        self.scale_dict = None
        self.symbolic_scale_dict = None

        # 2D array with number of models rejected from completions of jth
        # ordered operator because a subset features in completions of ith
        # ordered operator
        self.filter_data = np.zeros([len(self.data), len(self.data)])

    @property
    def is_democratic(self):
        return self.philosophy == "democratic"

    @property
    def is_stringent(self):
        return self.philosophy == "stringent"

    @property
    def is_mass(self):
        return self.criterion == "mass"

    @property
    def is_dimension(self):
        return self.criterion == "dimension"

    def query(self, func):
        """Function that acts on each model"""
        return {k: [m for m in v if func(m)] for k, v in self.data.items()}

    def filter_by_query(self, func):
        """Alter internal data from results of query by side effect"""
        self.data = self.query(func)

    @classmethod
    def no_seesaws(cls, model):
        """Aux. query to remove seesaw fields"""
        no_seesaws = (
            lambda m: not m.contains_field("F,00,2,0,0")
            and not m.contains_field("F,00,0,0,0")
            and not m.contains_field("S,00,2,1,0")
        )
        return no_seesaws(model)

    # model number is the product of primes representing the model
    def democratic_model_number(self, model):
        """Product of primes representing the fields in the model"""
        prod = 1
        for qn in model.quantum_numbers:
            prod *= self.exotic_prime_dict[qn]
        return prod

    def stringent_model_number(self, model):
        """Product of primes representing the terms in the Lagrangian of the model"""
        prod = 1
        for term in model.head["terms"]:
            prod *= self.term_prime_dict[term]
        return prod

    def model_number(self, model):
        """General dispatch for model number"""
        if self.is_democratic:
            return self.democratic_model_number(model)
        return self.stringent_model_number(model)

    def force(self):
        """Upgrade internal data from LazyCompletion objects to Completion objects. May
        take a while to run, probably filter the database down a little before
        running this.

        """
        if self.is_forced:
            return

        self.is_forced = True
        for k, v in self.data.items():
            self.data = {k: [m.force() for m in v]}

    def democratic_remove_equivalent_models(self):
        """Removes duplicate models only by field content"""
        from neutrinomass.utils.functions import remove_equivalent_nopop

        def eq(x, y):
            return self.democratic_model_number(x) == self.democratic_model_number(y)

        for k, v in self.data.items():
            self.data[k] = remove_equivalent_nopop(v, eq_func=eq)

    def stringent_remove_equivalent_models(self):
        """Removes duplicate models by interaction terms in the Lagrangian"""
        from neutrinomass.utils.functions import remove_equivalent_nopop

        def eq(x, y):
            return self.stringent_model_number(x) == self.stringent_model_number(y)

        for k, v in self.data.items():
            self.data[k] = remove_equivalent_nopop(v, eq_func=eq)

    def remove_equivalent_models(self):
        """General dispatch function for removing equivalent models"""
        if self.is_democratic:
            return self.democratic_remove_equivalent_models()
        return self.stringent_remove_equivalent_models()

    def filter_model_by_mass(self, op: str, model, keep_filter_data=False):
        """Remove all completions with the same or a subset of the particle content of
        an upstream model by the neutrino-mass criterion, i.e. only keep
        leading-order contributions to the neutrino mass.

        Dispatch on filtering philosophy handled by call to `model_number`.

        """
        op_scale = self.scale_dict[op]
        sieve = self.model_number(model)
        for k, v in self.data.items():
            if self.scale_dict[k] >= op_scale:
                continue

            new_v = []
            for test_model in v:
                if self.model_number(test_model) % sieve != 0:
                    new_v.append(test_model)
                    continue

                if keep_filter_data:
                    ordered_op_label_list = list(self.data)
                    sieve_op_pos = ordered_op_label_list.index(op)
                    other_op_pos = ordered_op_label_list.index(k)
                    self.filter_data[sieve_op_pos][other_op_pos] += 1 / len(v)

            self.data[k] = new_v

    def filter_one_loop_weinberg(self, democratic_nums):
        one_loop_scale = 605520000000.0 / (16 * math.pi ** 2)
        for k, v in self.data.items():
            if self.scale_dict[k] >= one_loop_scale:
                continue

            new_v = []
            for model in v:
                for n in democratic_nums:
                    if self.democratic_model_number(model) % n == 0:
                        keep = False
                        break
                else:
                    keep = True

                if keep:
                    new_v.append(model)

            self.data[k] = new_v

    def filter_models_by_mass(self, op: str):
        for model in self.data[op]:
            self.filter_model_by_mass(op, model)

    def filter_by_mass(self):
        for op in self.data:
            self.filter_models_by_mass(op)

    def filter_model_by_dimension(self, op: str, model):
        """Remove all completions with the same or a subset of the particle content of
        an upstream model by the dimension criterion, i.e. only keep models that
        would imply lower dimensional operators.

        Dispatch on filtering philosophy handled by call to `model_number`.

        """
        from neutrinomass.completions import EFF_OPERATORS
        from neutrinomass.completions import DERIV_EFF_OPERATORS

        ops = {**EFF_OPERATORS, **DERIV_EFF_OPERATORS}

        op_dim = ops[op].mass_dimension
        sieve = self.model_number(model)
        for k, v in self.data.items():
            if ops[k].mass_dimension <= op_dim:
                continue

            new_v = []
            for test_model in v:
                if self.model_number(test_model) % sieve != 0:
                    new_v.append(test_model)
                    continue

                ordered_op_label_list = list(self.data)
                sieve_op_pos = ordered_op_label_list.index(op)
                other_op_pos = ordered_op_label_list.index(k)
                self.filter_data[sieve_op_pos][other_op_pos] += 1 / len(v)

            self.data[k] = new_v

    def filter_models_by_dimension(self, op: str):
        for model in self.data[op]:
            self.filter_model_by_dimension(op, model)

    def filter_by_dimension(self):
        for op in self.data:
            self.filter_models_by_dimension(op)

    def filter(self):
        """Filter dispatch on filtering criterion"""
        if not self.is_ordered:
            self.order()

        if self.is_mass:
            self.filter_by_mass()
        else:
            self.filter_by_dimension()

    def fill_scale_dict(self):
        from neutrinomass.database import neutrino_mass_estimate
        from neutrinomass.database import numerical_np_scale_estimate
        from neutrinomass.completions import EFF_OPERATORS
        from neutrinomass.completions import DERIV_EFF_OPERATORS

        scale_dict, symbolic_scale_dict = {}, {}
        ops = {**EFF_OPERATORS, **DERIV_EFF_OPERATORS}
        for k, v in ops.items():
            estimates = neutrino_mass_estimate(v)
            max_scale = 0
            max_symbolic_scale = None
            for symbolic_expr in estimates:
                scale = numerical_np_scale_estimate(symbolic_expr)
                if scale > max_scale:
                    max_scale = scale
                    max_symbolic_scale = symbolic_expr

            assert max_scale > 0
            assert max_symbolic_scale is not None

            scale_dict[k] = max_scale
            symbolic_scale_dict[k] = max_symbolic_scale

        self.scale_dict = scale_dict
        self.symbolic_scale_dict = symbolic_scale_dict

    def order_by_mass(self):
        """Provides `scale_dict` and orders the data dictionary by neutrino mass scale
        prediction.

        Currently does the same as fill_scale_dict but also orders data

        """
        from neutrinomass.database import neutrino_mass_estimate
        from neutrinomass.database import numerical_np_scale_estimate
        from neutrinomass.completions import EFF_OPERATORS
        from neutrinomass.completions import DERIV_EFF_OPERATORS

        def scale_pred(op):
            return max(
                numerical_np_scale_estimate(i) for i in neutrino_mass_estimate(op)
            )

        ops = {**EFF_OPERATORS, **DERIV_EFF_OPERATORS}
        scales = {k: scale_pred(v) for k, v in ops.items() if k in self.data}
        mv_ordered = dict(reversed(sorted(scales.items(), key=lambda x: x[1])))

        self.scale_dict = {k: v for k, v in mv_ordered.items()}
        self.data = {k: self.data[k] for k, v in mv_ordered.items()}
        self.is_ordered = True

    def order_by_dimension(self):
        """Orders the data dictionary by operator dimension"""
        from neutrinomass.completions import EFF_OPERATORS
        from neutrinomass.completions import DERIV_EFF_OPERATORS

        ops = {**EFF_OPERATORS, **DERIV_EFF_OPERATORS}
        sorted_data = sorted(
            self.data.items(), key=lambda kv: ops[kv[0]].mass_dimension
        )

        self.data = dict(sorted_data)
        self.is_ordered = True

    def order(self):
        """General dispatch on order by filtering criterion"""
        if self.is_mass:
            self.order_by_mass()
        else:
            self.order_by_dimension()

    def process(self, filter_seesaws=False):
        """All common preprocessing steps in one function"""
        if filter_seesaws:
            # print("Removing seesaw fields...")
            self.filter_by_query(ModelDatabase.no_seesaws)
        # print("Removing equivalent models...")
        self.remove_equivalent_models()
        # print("Ordering by neutrino-mass scales...")
        self.order()
        # print("Filtering democratically by neutrino-mass estimate...")
        self.filter()


DATA = pickle.load(open(os.path.join(os.path.dirname(__file__), "democratic.p"), "rb"))
EXOTICS = pickle.load(open(os.path.join(os.path.dirname(__file__), "exotics.p"), "rb"))
TERMS = pickle.load(open(os.path.join(os.path.dirname(__file__), "terms.p"), "rb"))
MVDF = ModelDataFrame.new(data=DATA, exotics=EXOTICS, terms=TERMS)
