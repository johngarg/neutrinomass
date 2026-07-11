#!/usr/bin/env python3

import json

import networkx as nx
import pytest

from neutrinomass.completions.completions import (
    are_equivalent_completions,
    base_exotic_label,
    operator_completions,
)
from neutrinomass.completions.core import Completion, DerivativeRoute, LorentzProjection
from neutrinomass.completions.fingerprints import completion_fingerprint
from neutrinomass.completions.operators import EFF_OPERATORS, DERIV_EFF_OPERATORS
from neutrinomass.database import dumps_completion as public_dumps_completion
from neutrinomass.database.serialization import (
    completion_from_record,
    completion_to_record,
    dumps_completion,
    loads_completion,
    operator_from_data,
    operator_to_data,
    read_completion_jsonl,
    write_completion_jsonl,
)


def completion_with_route():
    completion = next(operator_completions(EFF_OPERATORS["1"]))
    exotic = next(iter(completion.exotics))
    edge = next(
        edge
        for edge, particle in nx.get_edge_attributes(
            completion.graph, "particle"
        ).items()
        if base_exotic_label(particle) == base_exotic_label(exotic.label)
    )
    route = DerivativeRoute(
        edge=edge,
        numerator_field=exotic.label,
        numerator_lorentz=exotic.dynkin[:2],
        differentiated_field="L",
        differentiated_lorentz="10",
    )
    return Completion(
        operator=completion.operator,
        partition=completion.partition,
        graph=completion.graph,
        exotics=completion.exotics,
        terms=completion.terms,
        topology=completion.topology,
        canonical_topology=completion.canonical_topology,
        derivative_routes=(route,),
    )


def test_completion_json_round_trip_is_exact_and_non_executable():
    completion = completion_with_route()
    payload = dumps_completion(completion)
    restored = loads_completion(payload)

    assert public_dumps_completion(completion) == payload
    assert "Completion(" not in payload
    assert completion_fingerprint(restored) == completion_fingerprint(completion)
    assert are_equivalent_completions(restored, completion)
    assert dumps_completion(restored) == payload
    assert restored.derivative_routes == completion.derivative_routes
    assert nx.to_dict_of_dicts(restored.graph) == nx.to_dict_of_dicts(
        completion.graph
    )


def test_derivative_operator_round_trip_preserves_stripped_metadata():
    operator = DERIV_EFF_OPERATORS["D20"].operator
    restored = operator_from_data(operator_to_data(operator))

    assert restored == operator
    assert [field.stripped for field in restored.indexed_fields] == [
        field.stripped for field in operator.indexed_fields
    ]


def test_lorentz_projection_round_trip_preserves_basis_metadata():
    completion = completion_with_route()
    completion.lorentz_projection = LorentzProjection(
        basis_labels=("u:0-1,2-3|d:0-1", "u:0-2,1-3|d:0-1"),
        coordinates=("1", "-1"),
        derivative_field="DH",
        ibp_relation="D(H1) H2 + H1 D(H2) = 0 modulo a total derivative",
        eom_relation="derivatives on external fermions are removed by their EOM",
    )

    restored = loads_completion(dumps_completion(completion))

    assert restored.lorentz_projection == completion.lorentz_projection
    assert completion_fingerprint(restored) == completion_fingerprint(completion)


def test_completion_jsonl_round_trip(tmp_path):
    completions = [
        completion_with_route(),
        next(operator_completions(EFF_OPERATORS["2"])),
    ]
    path = tmp_path / "completions.jsonl"

    write_completion_jsonl(path, completions)
    restored = read_completion_jsonl(path)

    assert [completion_fingerprint(item) for item in restored] == [
        completion_fingerprint(item) for item in completions
    ]


def test_completion_schema_rejects_unknown_versions():
    record = completion_to_record(completion_with_route())
    record["version"] += 1

    with pytest.raises(ValueError, match="Unsupported completion schema"):
        completion_from_record(json.loads(json.dumps(record)))


def test_completion_schema_rejects_executable_charge_expressions():
    record = completion_to_record(completion_with_route())
    record["operator"]["tensors"][0]["charges"]["y"] = (
        "__import__('os').system('echo unsafe')"
    )

    with pytest.raises(ValueError, match="Invalid rational charge"):
        completion_from_record(record)


def test_completion_schema_rejects_executable_projection_coordinates():
    completion = completion_with_route()
    completion.lorentz_projection = LorentzProjection(
        basis_labels=("basis-0",),
        coordinates=("1",),
        derivative_field="DH",
        ibp_relation="test relation",
        eom_relation="test relation",
    )
    record = completion_to_record(completion)
    record["lorentz_projection"]["coordinates"][0] = "__import__('os').system('x')"

    with pytest.raises(ValueError, match="Invalid Lorentz-projection coordinate"):
        completion_from_record(record)
