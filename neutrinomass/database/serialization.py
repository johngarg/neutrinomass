#!/usr/bin/env python3

"""Versioned, non-executable JSON serialisation for completion records."""

from pathlib import Path
import json
import re

import networkx as nx
from sympy import Rational

from neutrinomass.completions.core import (
    Completion,
    ComplexScalar,
    DerivativeRoute,
    EffectiveOperator,
    FieldType,
    MajoranaFermion,
    RealScalar,
    VectorLikeDiracFermion,
)
from neutrinomass.completions.topologies import Leaf
from neutrinomass.tensormethod.core import IndexedField, Operator, delta, eps


SCHEMA_NAME = "neutrinomass.completion"
SCHEMA_VERSION = 1


_RATIONAL_PATTERN = re.compile(r"[+-]?\d+(?:/[1-9]\d*)?\Z")


def _decode_charge(value):
    """Decode the schema's deliberately restricted integer/rational grammar."""

    if not isinstance(value, str) or _RATIONAL_PATTERN.fullmatch(value) is None:
        raise ValueError(f"Invalid rational charge {value!r}")
    return Rational(value)


def _encode_stripped(stripped):
    if stripped is None:
        return None
    encoded = dict(stripped)
    if encoded.get("charges") is not None:
        encoded["charges"] = {
            name: str(value) for name, value in encoded["charges"].items()
        }
    return encoded


def _decode_stripped(stripped):
    if stripped is None:
        return None
    decoded = dict(stripped)
    if decoded.get("charges") is not None:
        decoded["charges"] = {
            name: _decode_charge(value)
            for name, value in decoded["charges"].items()
        }
    return decoded


def tensor_to_data(tensor):
    """Convert an indexed field or invariant tensor to JSON-compatible data."""

    if isinstance(tensor, IndexedField):
        field_class = (
            tensor.__class__.__name__
            if isinstance(tensor, FieldType)
            else "IndexedField"
        )
        return {
            "kind": "field",
            "field_class": field_class,
            "label": tensor.label,
            "indices": [str(index) for index in tensor.indices],
            "charges": {
                name: str(value) for name, value in tensor.charges.items()
            },
            "is_conj": tensor.is_conj,
            "symmetry": tensor.symmetry,
            "comm": tensor.comm,
            "latex": tensor.latex,
            "nf": tensor.nf,
            "derivs": tensor.derivs,
            "stripped": _encode_stripped(tensor.stripped),
            "is_unbarred": getattr(tensor, "is_unbarred", None),
        }

    label = str(tensor).partition("(")[0]
    if label in {"metric", "Eps"}:
        kind = "epsilon"
    elif label == "KD":
        kind = "delta"
    else:
        raise ValueError(f"Unrecognised tensor {tensor}")
    return {
        "kind": kind,
        "indices": [str(index) for index in tensor.indices],
    }


def tensor_from_data(data):
    """Reconstruct a tensor from versioned JSON data."""

    kind = data["kind"]
    indices = " ".join(data["indices"])
    if kind == "epsilon":
        return eps(indices)
    if kind == "delta":
        return delta(indices)
    if kind != "field":
        raise ValueError(f"Unknown tensor kind {kind}")

    charges = {
        name: _decode_charge(value) for name, value in data["charges"].items()
    }
    field_class = data["field_class"]
    common = {
        "label": data["label"],
        "indices": indices,
        "charges": charges,
        "latex": data["latex"],
        "is_conj": data["is_conj"],
        "symmetry": data["symmetry"],
    }
    if field_class == "ComplexScalar":
        return ComplexScalar(**common)
    if field_class == "RealScalar":
        return RealScalar(**common)
    if field_class == "MajoranaFermion":
        return MajoranaFermion(**common)
    if field_class == "VectorLikeDiracFermion":
        return VectorLikeDiracFermion(
            **common, is_unbarred=data["is_unbarred"]
        )
    if field_class != "IndexedField":
        raise ValueError(f"Unknown field class {field_class}")
    return IndexedField(
        **common,
        comm=data["comm"],
        nf=data["nf"],
        derivs=data["derivs"],
        stripped=_decode_stripped(data["stripped"]),
    )


def operator_to_data(operator):
    return [tensor_to_data(tensor) for tensor in operator.tensors]


def operator_from_data(data):
    return Operator(*(tensor_from_data(tensor) for tensor in data))


def partition_to_data(partition):
    if isinstance(partition, Leaf):
        return {
            "kind": "leaf",
            "field": tensor_to_data(partition.field),
            "node": partition.node,
        }
    return {
        "kind": "branch",
        "children": [partition_to_data(child) for child in partition],
    }


def partition_from_data(data):
    if data["kind"] == "leaf":
        return Leaf(tensor_from_data(data["field"]), data["node"])
    if data["kind"] != "branch":
        raise ValueError(f"Unknown partition kind {data['kind']}")
    return tuple(partition_from_data(child) for child in data["children"])


def graph_to_data(graph):
    edges = sorted(
        (min(source, target), max(source, target), dict(attributes))
        for source, target, attributes in graph.edges(data=True)
    )
    return {
        "nodes": [
            {"id": node, "attributes": dict(attributes)}
            for node, attributes in sorted(graph.nodes(data=True), key=lambda item: item[0])
        ],
        "edges": [
            {
                "source": source,
                "target": target,
                "attributes": dict(attributes),
            }
            for source, target, attributes in edges
        ],
    }


def graph_from_data(data):
    graph = nx.Graph()
    for node in data["nodes"]:
        graph.add_node(node["id"], **node["attributes"])
    for edge in data["edges"]:
        graph.add_edge(
            edge["source"], edge["target"], **edge["attributes"]
        )
    return graph


def route_to_data(route):
    return {
        "edge": list(route.edge),
        "numerator_field": route.numerator_field,
        "numerator_lorentz": route.numerator_lorentz,
        "differentiated_field": route.differentiated_field,
        "differentiated_lorentz": route.differentiated_lorentz,
    }


def route_from_data(data):
    return DerivativeRoute(
        edge=tuple(data["edge"]),
        numerator_field=data["numerator_field"],
        numerator_lorentz=data["numerator_lorentz"],
        differentiated_field=data["differentiated_field"],
        differentiated_lorentz=data["differentiated_lorentz"],
    )


def completion_to_record(completion):
    return {
        "schema": SCHEMA_NAME,
        "version": SCHEMA_VERSION,
        "operator": {
            "name": completion.operator.name,
            "tensors": operator_to_data(completion.operator.operator),
        },
        "partition": partition_to_data(completion.partition),
        "graph": graph_to_data(completion.graph),
        "exotics": [
            tensor_to_data(field)
            for field in sorted(completion.exotics, key=lambda field: field.label)
        ],
        "terms": [operator_to_data(term) for term in completion.terms],
        "topology": completion.topology,
        "canonical_topology": completion.canonical_topology,
        "derivative_routes": [
            route_to_data(route) for route in completion.derivative_routes
        ],
        "legacy_derivative_edges": [
            list(edge)
            for edge in completion.derivative_edges
            if not completion.derivative_routes
        ],
    }


def completion_from_record(record):
    if record.get("schema") != SCHEMA_NAME:
        raise ValueError("Not a neutrinomass completion record")
    if record.get("version") != SCHEMA_VERSION:
        raise ValueError(f"Unsupported completion schema {record.get('version')}")

    operator = EffectiveOperator(
        record["operator"]["name"],
        operator_from_data(record["operator"]["tensors"]),
    )
    routes = tuple(route_from_data(route) for route in record["derivative_routes"])
    return Completion(
        operator=operator,
        partition=partition_from_data(record["partition"]),
        graph=graph_from_data(record["graph"]),
        exotics={tensor_from_data(field) for field in record["exotics"]},
        terms=[operator_from_data(term) for term in record["terms"]],
        topology=record["topology"],
        canonical_topology=record["canonical_topology"],
        derivative_edges=tuple(
            tuple(edge) for edge in record["legacy_derivative_edges"]
        ),
        derivative_routes=routes,
    )


def dumps_completion(completion):
    return json.dumps(
        completion_to_record(completion),
        sort_keys=True,
        separators=(",", ":"),
    )


def loads_completion(payload):
    return completion_from_record(json.loads(payload))


def write_completion_jsonl(path, completions):
    path = Path(path)
    with path.open("w", encoding="utf-8") as stream:
        for completion in completions:
            stream.write(dumps_completion(completion) + "\n")


def read_completion_jsonl(path):
    path = Path(path)
    with path.open("r", encoding="utf-8") as stream:
        return [loads_completion(line) for line in stream if line.strip()]
