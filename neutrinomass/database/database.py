#!/usr/bin/env python3

from neutrinomass.tensormethod import D, L, Q, H, eb, ub, db
from neutrinomass.completions.core import Completion, cons_completion_field
from networkx import Graph

ExoticField = cons_completion_field


def Partition(*args):
    return tuple(args)


def read_completion(comp: str):
    pass


def read_completions(filename: str):
    """Do this as a context manager?"""
    pass
