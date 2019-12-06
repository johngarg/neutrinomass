#!/usr/bin/env python3

"""Functions to read and parse topologies generated with Mathematica."""

import os
from glob import glob
import matplotlib.pyplot as plt
import networkx as nx

TOPOLOGY_PATH = "/Users/johngargalionis/Dropbox/PhD/mv/mv/mv/completions/topology_data/"
PARTITIONS = os.path.join(TOPOLOGY_PATH, "partitions")
DIAGRAMS = os.path.join(TOPOLOGY_PATH, "diagrams")
GRAPHS = os.path.join(TOPOLOGY_PATH, "graphs")


def read_topology_file(data_path) -> str:
    """Reads the topology and returns the contents of the data file as a string."""
    with open(data_path, "r") as f:
        data_string = f.read()

    return data_string


def eval_partition(partition: str):
    S = lambda x: ("S", x)
    F = lambda x: ("F", x)

    def List(*args):
        return args

    structure = eval(partition)

    # Take first element to simplify output but ensure not losing any info
    assert len(structure) == 1
    return structure[0]


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
        # datum["diagram"] = img
        out.append(topology)

    return out
