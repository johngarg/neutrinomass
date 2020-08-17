#!/usr/bin/env python3

from neutrinomass.tensormethod import D, L, Q, H, eb, ub, db, eps, delta
from neutrinomass.tensormethod.core import IndexedField, Field
from neutrinomass.completions.topologies import Leaf
from neutrinomass.completions.core import (
    EffectiveOperator,
    Completion,
    cons_completion_field,
)
from networkx import Graph
from collections import defaultdict
import time
from typing import Dict, List
import re

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
        return self.head["quantum_numbers"]

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


class ModelDatabase:
    def __init__(self, path: str):
        self.path = path

        import os
        from glob import glob

        print("Initialising database...")
        filenames = glob(os.path.join(self.path, "*.dat"))
        mvdb_data = [dict(read_completions(f)) for f in filenames]
        mvdb = {k: v for d in mvdb_data for k, v in d.items()}
        self.data = mvdb

    def query(self, func):
        """Function that acts on each model"""
        return {k: [m for m in v if func(m)] for k, v in self.data.items()}

    def filter_seesaw(self):
        no_seesaws = (
            lambda m: not m.contains_field("F,00,2,0,0")
            and not m.contains_field("F,00,0,0,0")
            and not m.contains_field("S,00,2,1,0")
        )
        return self.query(no_seesaws)
