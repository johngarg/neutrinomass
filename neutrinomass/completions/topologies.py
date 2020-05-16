#!/usr/bin/env python3

"""Functions to read and parse topologies generated with Mathematica."""

import os
from glob import glob
import matplotlib.pyplot as plt
import networkx as nx
from typing import NamedTuple
from neutrinomass.tensormethod.core import IndexedField

# PATH_TO_MV = "/Users/johngargalionis/Dropbox/PhD/mv/"
TOPOLOGY_PATH = os.path.join(os.path.dirname(__file__), "topology_data")
# INTERNAL_PATH = "neutrinomass/neutrinomass/completions/topology_data/"
# TOPOLOGY_PATH = PATH_TO_MV + INTERNAL_PATH
PARTITIONS = os.path.join(TOPOLOGY_PATH, "partitions")
DIAGRAMS = os.path.join(TOPOLOGY_PATH, "diagrams")
GRAPHS = os.path.join(TOPOLOGY_PATH, "graphs")


class Leaf(NamedTuple):
    field: IndexedField
    node: int


def read_topology_file(data_path) -> str:
    """Reads the topology and returns the contents of the data file as a string."""
    with open(data_path, "r") as f:
        data_string = f.read()

    return data_string


def eval_partition(partition: str):
    S = lambda x: Leaf("S", x)
    F = lambda x: Leaf("F", x)

    def List(*args):
        return args

    structure = eval(partition)

    # Take first element to simplify output but ensure not losing any info
    return structure


def eval_graph(graph: str):
    G = nx.Graph()
    for edge in graph.splitlines():
        i, j = eval(edge)
        G.add_edge(i, j)

    return G


def get_topology_data(n_scalars, n_fermions):
    """Returns a list of dictionaries with data from topology data files.

    [{"partition": parition_string, "graph": graph_string, "img": image}]

    """
    partition_files = sorted(glob(PARTITIONS + f"/{n_scalars}s{n_fermions}f_*"))
    diagram_files = sorted(glob(DIAGRAMS + f"/{n_scalars}s{n_fermions}f_*"))
    graph_files = sorted(glob(GRAPHS + f"/{n_scalars}s{n_fermions}f_*"))

    if not partition_files:
        raise Exception("Topologies not found, please generate them again.")

    out = []
    for p, d, g in zip(partition_files, diagram_files, graph_files):
        topology = {}
        partition_string = eval_partition(read_topology_file(p))
        # img = plt.imread(d)
        graph_string = eval_graph(read_topology_file(g))

        topology["partition"] = partition_string
        topology["graph"] = graph_string
        # topology["diagram"] = img
        topology["partition_file"] = p
        out.append(topology)

    return out
