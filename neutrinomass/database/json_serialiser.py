#!/usr/bin/env python3

from neutrinomass.completions import *
from neutrinomass.completions.topologies import Leaf
from neutrinomass.completions.core import FieldType
from neutrinomass.tensormethod import Field

import networkx as nx
import json


def export_tensor(tensor):
    if isinstance(tensor, FieldType):
        indices = {k: [str(i) for i in v] for (k, v) in tensor.indices_by_type.items()}
        charges = {k: str(v) for k, v in tensor.charges.items()}
        kwargs = {
            "Label": tensor.label,
            "Indices": indices,
            "Symmetry": tensor.symmetry,
            "Charges": charges,
            "Nf": tensor.nf,
            "LorentzDynkin": tensor.lorentz_irrep,
            "IsospinDynkin": tensor.isospin_irrep[0],
            "ColourDynkin": tensor.colour_irrep,
            "Comm": tensor.comm,
            "Latex": tensor.latex,
            "IsConj": tensor.is_conj,
            "Exotic": True,
            "Field": True,
            "Deriv": False,
        }

        is_unbarred = None
        if hasattr(tensor, "is_unbarred"):
            is_unbarred = tensor.is_unbarred
        kwargs["IsUnbarred"] = is_unbarred
        return kwargs

    if isinstance(tensor, Field):
        # SM field
        indices = {k: [str(i) for i in v] for (k, v) in tensor.indices_by_type.items()}
        if tensor.derivs:
            assert tensor.derivs == 1
            a, b = tensor.lorentz_irrep
            dynkin_label = str(a) + str(b)

            field = tensor.strip_derivs()
            label = field.label
            return {
                "Label": label,
                "Indices": indices,
                "IsConj": tensor.is_conj,
                "DerivDynkin": dynkin_label,
                "Exotic": False,
                "Field": True,
                "Deriv": True,
            }

        label = tensor.field.label
        return {
            "Label": label,
            "Indices": indices,
            "IsConj": tensor.is_conj,
            "Field": True,
            "Exotic": False,
            "Deriv": False,
        }

    indices = [str(i) for i in tensor.indices]
    if str(tensor).startswith("metric") or str(tensor).startswith("Eps"):
        return {"Field": False, "label": "Eps", "Indices": indices}

    if str(tensor).startswith("KD"):
        return {"Field": False, "label": "Delta", "Indices": indices}

    raise ValueError(f"Unrecognised tensor: {tensor}")


def export_eff_operator(effop):
    return {"Name": effop.name, "Operator": export_operator(effop.operator)}


def export_partition(expr):
    if isinstance(expr, Leaf):
        return {"Field": str(expr.field), "Node": expr.node}

    return tuple([export_partition(i) for i in expr])


def export_graph(g: nx.Graph):
    data = nx.readwrite.json_graph.node_link_data(g)
    for link in data["links"]:
        link["particle"] = export_tensor(link["particle"])

    return data


def export_operator(op):
    tensors = op.tensors
    return [export_tensor(t) for t in op.tensors]


def export_terms(terms):
    return [export_operator(op) for op in terms]


def export_exotics(expr: set):
    return [export_tensor(f) for f in expr]


def export_completion(comp):
    """Function to export completions to json"""
    return {
        "EffectiveOperator": export_eff_operator(comp.operator),
        "Partition": export_partition(comp.partition),
        "Graph": export_graph(comp.graph),
        "Exotics": export_exotics(comp.exotics),
        "Terms": export_terms(comp.terms),
    }
