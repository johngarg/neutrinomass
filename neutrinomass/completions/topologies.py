#!/usr/bin/env python3

"""Functions to read and parse topologies generated with Mathematica."""

import os
from glob import glob
from pathlib import Path
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


# The 5s2f panels in Fig. 10 of arXiv:2009.13537 were assembled with the
# topology files in lexicographic filename order and then labelled 1--24.
# Keep that published numbering as the user-facing label while retaining the
# generated filename as the canonical topology identifier.
PAPER_TOPOLOGY_PANEL_ORDER = {
    "5s2f": (
        1,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        2,
        20,
        21,
        22,
        23,
        24,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
    )
}


class Leaf(NamedTuple):
    field: IndexedField
    node: int


def split_topology_name(topology):
    """Return the field-content prefix and integer suffix of a topology name."""

    prefix, number = Path(topology).stem.rsplit("_", 1)
    return prefix, int(number)


def paper_topology_name(topology):
    """Translate a generated topology identifier to its published label."""

    prefix, number = split_topology_name(topology)
    panel_order = PAPER_TOPOLOGY_PANEL_ORDER.get(prefix)
    if panel_order is None:
        return f"{prefix}_{number}"

    try:
        paper_number = panel_order.index(number) + 1
    except ValueError as error:
        raise ValueError(f"Unknown canonical topology: {prefix}_{number}") from error

    return f"{prefix}_{paper_number}"


def canonical_topology_name(topology):
    """Translate a published topology label to its generated identifier."""

    prefix, paper_number = split_topology_name(topology)
    panel_order = PAPER_TOPOLOGY_PANEL_ORDER.get(prefix)
    if panel_order is None:
        return f"{prefix}_{paper_number}"
    if not 1 <= paper_number <= len(panel_order):
        raise ValueError(f"Unknown published topology: {prefix}_{paper_number}")

    return f"{prefix}_{panel_order[paper_number - 1]}"


def paper_topology_sort_key(topology):
    """Sort topology paths by the numerical label used in the paper."""

    return split_topology_name(paper_topology_name(topology))


def topology_files_by_stem(directory, pattern):
    """Return topology files keyed by their generated topology identifier."""

    return {Path(path).stem: path for path in glob(os.path.join(directory, pattern))}


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
    pattern = f"{n_scalars}s{n_fermions}f_*"
    partition_files = topology_files_by_stem(PARTITIONS, pattern)
    diagram_files = topology_files_by_stem(DIAGRAMS, pattern)
    graph_files = topology_files_by_stem(GRAPHS, pattern)

    if not partition_files:
        raise Exception("Topologies not found, please generate them again.")

    file_stems = set(partition_files)
    if file_stems != set(diagram_files) or file_stems != set(graph_files):
        raise Exception("Topology partitions, diagrams and graphs do not match.")

    out = []
    for canonical_name in sorted(file_stems, key=paper_topology_sort_key):
        p = partition_files[canonical_name]
        d = diagram_files[canonical_name]
        g = graph_files[canonical_name]
        topology = {}
        partition_string = eval_partition(read_topology_file(p))
        # img = plt.imread(d)
        graph_string = eval_graph(read_topology_file(g))

        topology["partition"] = partition_string
        topology["graph"] = graph_string
        # topology["diagram"] = img
        topology["topology"] = paper_topology_name(canonical_name)
        topology["canonical_topology"] = canonical_name
        topology["partition_file"] = p
        topology["diagram_file"] = d
        topology["graph_file"] = g
        out.append(topology)

    return out
