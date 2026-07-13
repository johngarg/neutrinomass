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
    LorentzProjection,
    MajoranaFermion,
    MultiDerivativeProjection,
    PropagatorContribution,
    RealScalar,
    VectorLikeDiracFermion,
)
from neutrinomass.completions.topologies import Leaf
from neutrinomass.tensormethod.core import IndexedField, Operator, delta, eps


SCHEMA_NAME = "neutrinomass.completion"
SCHEMA_VERSION = 2


_RATIONAL_PATTERN = re.compile(r"[+-]?\d+(?:/[1-9]\d*)?\Z")


class CompletionJSONLError(ValueError):
    """Report a malformed completion record at its physical JSONL line."""

    def __init__(self, path, line_number, error):
        self.path = Path(path)
        self.line_number = line_number
        self.error = error
        super().__init__(f"{self.path}:{line_number}: {error}")


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


def contribution_to_data(contribution):
    return {
        "edge": list(contribution.edge),
        "particle": contribution.particle,
        "numerator_kind": contribution.numerator_kind,
        "denominator_order": contribution.denominator_order,
        "numerator_lorentz": contribution.numerator_lorentz,
        "differentiated_field": contribution.differentiated_field,
        "differentiated_lorentz": contribution.differentiated_lorentz,
        "cut_side": list(contribution.cut_side),
    }


def contribution_from_data(data):
    return PropagatorContribution(
        edge=tuple(data["edge"]),
        particle=data["particle"],
        numerator_kind=data["numerator_kind"],
        denominator_order=data["denominator_order"],
        numerator_lorentz=data["numerator_lorentz"],
        differentiated_field=data["differentiated_field"],
        differentiated_lorentz=data["differentiated_lorentz"],
        cut_side=tuple(data["cut_side"]),
    )


def projection_to_data(projection):
    if projection is None:
        return None
    data = {
        "basis_labels": list(projection.basis_labels),
        "coordinates": list(projection.coordinates),
        "ibp_relation": projection.ibp_relation,
        "eom_relation": projection.eom_relation,
    }
    if isinstance(projection, MultiDerivativeProjection):
        data["kind"] = "multi_derivative"
        data["derivative_fields"] = list(projection.derivative_fields)
    else:
        data["kind"] = "single_derivative"
        data["derivative_field"] = projection.derivative_field
    return data


def projection_from_data(data):
    if data is None:
        return None
    coordinates = tuple(data["coordinates"])
    if any(
        not isinstance(value, str)
        or _RATIONAL_PATTERN.fullmatch(value) is None
        for value in coordinates
    ):
        raise ValueError("Invalid Lorentz-projection coordinate")
    common = {
        "basis_labels": tuple(data["basis_labels"]),
        "coordinates": coordinates,
        "ibp_relation": data["ibp_relation"],
        "eom_relation": data["eom_relation"],
    }
    if data.get("kind", "single_derivative") == "multi_derivative":
        return MultiDerivativeProjection(
            derivative_fields=tuple(data["derivative_fields"]),
            **common,
        )
    return LorentzProjection(
        derivative_field=data["derivative_field"], **common
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
        "momentum_contributions": [
            contribution_to_data(contribution)
            for contribution in completion.momentum_contributions
        ],
        "lorentz_projection": projection_to_data(completion.lorentz_projection),
        "legacy_derivative_edges": [
            list(edge)
            for edge in completion.derivative_edges
            if not completion.derivative_routes
        ],
    }


def completion_from_record(record):
    if record.get("schema") != SCHEMA_NAME:
        raise ValueError("Not a neutrinomass completion record")
    version = record.get("version")
    if version not in {1, SCHEMA_VERSION}:
        raise ValueError(f"Unsupported completion schema {record.get('version')}")

    operator = EffectiveOperator(
        record["operator"]["name"],
        operator_from_data(record["operator"]["tensors"]),
    )
    if version == 1:
        routes = tuple(
            route_from_data(route) for route in record["derivative_routes"]
        )
        contributions = ()
    else:
        routes = ()
        contributions = tuple(
            contribution_from_data(contribution)
            for contribution in record["momentum_contributions"]
        )
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
        momentum_contributions=contributions,
        lorentz_projection=projection_from_data(record.get("lorentz_projection")),
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
    return list(iter_completion_jsonl(path))


def iter_completion_jsonl(path):
    """Yield completions from a JSONL artifact without retaining the full file."""

    path = Path(path)
    with path.open("r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, start=1):
            if not line.strip():
                continue
            try:
                yield loads_completion(line)
            except (
                AttributeError,
                IndexError,
                KeyError,
                TypeError,
                ValueError,
            ) as error:
                raise CompletionJSONLError(path, line_number, error) from error
