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

ExoticField = cons_completion_field


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

    # @property
    # def operator_dimension(self):
    #     return self.head["operator_dimension"]


def read_completions(filename: str):
    """Do this as a context manager?"""
    completions = defaultdict(list)
    with open(filename, "r") as f:
        line = f.readline()
        counter = 1
        while line:
            # start = time.time()
            comp = eval(line)
            completions[comp.operator_name].append(comp)
            line = f.readline()
            # print(f"Took {time.time() - start} seconds")
            counter += 1
            if counter == 1000:
                break

    return completions


def apply_naive_filter(
    db: Dict[int, List[LazyCompletion]]
) -> Dict[int, List[LazyCompletion]]:
    for mass_dim, comps in db.items():
        for comp in db:
            pass


def query(completions: List[LazyCompletion], func):
    for comp in completions:
        if func(comp):
            yield comp
