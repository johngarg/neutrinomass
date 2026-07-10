#!/usr/bin/env python3

"""Exact contraction-graph comparisons for completion interactions."""

from collections import Counter, OrderedDict

import networkx as nx

from neutrinomass.tensormethod.core import IndexedField


def field_label_parts(label: str):
    """Return a base label and its conjugate/Dirac-partner suffix bits."""

    suffix = ""
    while label and label[-1] in "†~":
        suffix = label[-1] + suffix
        label = label[:-1]
    return label, "†" in suffix, "~" in suffix


def _normalise_relabeling(relabeling):
    if isinstance(relabeling, str):
        return relabeling, False, False
    return relabeling


def field_relabeling(label: str, label_mapping):
    """Return the mapped label and the conjugate/Dirac orientation flips."""

    base_label, is_conjugate, is_dirac_partner = field_label_parts(label)
    target_label, conjugate_flip, dirac_flip = _normalise_relabeling(
        label_mapping.get(base_label, base_label)
    )
    target_label += "†" if is_conjugate ^ conjugate_flip else ""
    target_label += "~" if is_dirac_partner ^ dirac_flip else ""
    return target_label, conjugate_flip, dirac_flip


def remapped_field_label(label: str, label_mapping) -> str:
    """Relabel an exotic field, applying any supplied orientation flips."""

    return field_relabeling(label, label_mapping)[0]


def _index_key(index):
    raised_index = index if index.is_up else -index
    return index.index_type, str(raised_index)


def _field_signature(field, label_mapping):
    label = remapped_field_label(field.label, label_mapping)
    charges = tuple(
        sorted(
            (name, str(value)) for name, value in field.charges.items()
        )
    )
    return (
        "field",
        label,
        field.dynkin,
        charges,
        field.comm,
        field.derivs,
    )


def _invariant_signature(invariant):
    head = str(invariant).partition("(")[0]
    return "invariant", head, len(invariant.indices)


def _mapping_items(label_mapping):
    return tuple(
        (source_label, *_normalise_relabeling(relabeling))
        for source_label, relabeling in sorted(label_mapping.items())
    )


def _interaction_graph(term, mapping_items) -> nx.Graph:
    """Return the tensor-factor/index incidence graph for one interaction.

    Factor ordering and dummy-index names are deliberately absent from the graph.
    Index variance and invariant tensors remain explicit, so inequivalent gauge and
    Lorentz contractions do not collapse to the same representation.
    """

    cache_key = (id(term), mapping_items)
    cached = _GRAPH_CACHE.get(cache_key)
    if cached is not None and cached[0] is term:
        _GRAPH_CACHE.move_to_end(cache_key)
        return cached[1]

    label_mapping = {
        source_label: (target_label, conjugate_flip, dirac_flip)
        for source_label, target_label, conjugate_flip, dirac_flip in mapping_items
    }
    tensors = [tensor for tensor in term.tensors if hasattr(tensor, "indices")]
    index_counts = Counter(
        _index_key(index) for tensor in tensors for index in tensor.indices
    )
    index_nodes = {}
    graph = nx.Graph()

    for factor_number, tensor in enumerate(tensors):
        factor_node = ("factor", factor_number)
        if isinstance(tensor, IndexedField):
            signature = _field_signature(tensor, label_mapping)
        else:
            signature = _invariant_signature(tensor)
        graph.add_node(factor_node, signature=signature)

        for slot_number, index in enumerate(tensor.indices):
            index_key = _index_key(index)
            index_node = index_nodes.get(index_key)
            if index_node is None:
                index_node = ("index", len(index_nodes))
                index_nodes[index_key] = index_node
                graph.add_node(
                    index_node,
                    signature=(
                        "index",
                        index.index_type,
                        "free" if index_counts[index_key] == 1 else "contracted",
                    ),
                )

            port_node = ("port", factor_number, slot_number)
            graph.add_node(
                port_node,
                signature=("port", index.index_type, index.is_up),
            )
            graph.add_edge(factor_node, port_node)
            graph.add_edge(port_node, index_node)

    node_counts = Counter(
        (data["signature"], graph.degree[node])
        for node, data in graph.nodes(data=True)
    )
    graph.graph["coarse_key"] = tuple(sorted(node_counts.items(), key=repr))
    _GRAPH_CACHE[cache_key] = (term, graph)
    if len(_GRAPH_CACHE) > _GRAPH_CACHE_MAXSIZE:
        _GRAPH_CACHE.popitem(last=False)
    return graph


def interaction_graph(term, label_mapping=None) -> nx.Graph:
    """Return the tensor-factor/index incidence graph for one interaction.

    Factor ordering and dummy-index names are deliberately absent from the graph.
    Index variance and invariant tensors remain explicit, so inequivalent gauge and
    Lorentz contractions do not collapse to the same representation.
    """

    return _interaction_graph(term, _mapping_items(label_mapping or {}))


_GRAPH_CACHE_MAXSIZE = 50_000
_GRAPH_CACHE = OrderedDict()


_NODE_MATCH = nx.algorithms.isomorphism.categorical_node_match("signature", None)


def equivalent_interactions(term1, term2, label_mapping=None) -> bool:
    """Return whether two interactions have isomorphic contraction graphs."""

    graph1 = interaction_graph(term1, label_mapping)
    graph2 = interaction_graph(term2)
    if graph1.graph["coarse_key"] != graph2.graph["coarse_key"]:
        return False
    return nx.is_isomorphic(graph1, graph2, node_match=_NODE_MATCH)


def equivalent_lagrangians(terms1, terms2, label_mapping=None) -> bool:
    """Compare two interaction multisets, including repeated terms."""

    if len(terms1) != len(terms2):
        return False

    candidate_matches = nx.Graph()
    left_nodes = [("left", index) for index in range(len(terms1))]
    right_nodes = [("right", index) for index in range(len(terms2))]
    candidate_matches.add_nodes_from(left_nodes, bipartite=0)
    candidate_matches.add_nodes_from(right_nodes, bipartite=1)

    left_graphs = [interaction_graph(term, label_mapping) for term in terms1]
    right_graphs = [interaction_graph(term) for term in terms2]
    for left_number, left_graph in enumerate(left_graphs):
        for right_number, right_graph in enumerate(right_graphs):
            if left_graph.graph["coarse_key"] != right_graph.graph["coarse_key"]:
                continue
            if nx.is_isomorphic(left_graph, right_graph, node_match=_NODE_MATCH):
                candidate_matches.add_edge(
                    ("left", left_number), ("right", right_number)
                )

    matching = nx.algorithms.bipartite.maximum_matching(
        candidate_matches, top_nodes=left_nodes
    )
    return all(node in matching for node in left_nodes)
