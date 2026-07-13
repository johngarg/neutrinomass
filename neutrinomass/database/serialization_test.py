#!/usr/bin/env python3

import json

import networkx as nx
import pytest

from neutrinomass.completions.completions import (
    UnreducedSecondDerivativeBasis,
    are_equivalent_completions,
    base_exotic_label,
    expand_propagator_denominators,
    operator_completions,
)
from neutrinomass.completions.core import (
    Completion,
    DerivativeRoute,
    LorentzProjection,
    MultiDerivativeProjection,
    PropagatorContribution,
)
from neutrinomass.completions.fingerprints import completion_fingerprint
from neutrinomass.completions.operators import EFF_OPERATORS, DERIV_EFF_OPERATORS
from neutrinomass.database import dumps_completion as public_dumps_completion
from neutrinomass.database.serialization import (
    CompletionJSONLError,
    completion_from_record,
    completion_to_record,
    dumps_completion,
    iter_completion_jsonl,
    loads_completion,
    operator_from_data,
    operator_to_data,
    read_completion_jsonl,
    write_completion_jsonl,
    route_to_data,
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
    assert restored.momentum_contributions == completion.momentum_contributions
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


def test_multi_derivative_projection_round_trip_preserves_fields():
    completion = completion_with_route()
    completion.lorentz_projection = MultiDerivativeProjection(
        basis_labels=("u:0-1|d:0-1",),
        coordinates=("1",),
        derivative_fields=("DL", "DH"),
        ibp_relation="retain the named placement",
        eom_relation="No equation-of-motion reduction is applied",
    )

    restored = loads_completion(dumps_completion(completion))

    assert restored.lorentz_projection == completion.lorentz_projection
    assert restored.lorentz_projection.derivative_field == "DL,DH"
    assert completion_fingerprint(restored) == completion_fingerprint(completion)


def test_unreduced_d_squared_field_round_trip_preserves_derivative_count():
    basis = UnreducedSecondDerivativeBasis.from_operator(
        DERIV_EFF_OPERATORS["D15"]
    )
    operator = basis.placement_by_indices[(0, 0)].operator.operator

    restored = operator_from_data(operator_to_data(operator))
    boxed = next(field for field in restored.indexed_fields if field.derivs)

    assert boxed.derivs == 2
    assert boxed.strip_derivs().dynkin == boxed.dynkin


def test_multiple_propagator_contributions_round_trip():
    base = next(operator_completions(EFF_OPERATORS["2"]))
    completion = next(
        item
        for item in expand_propagator_denominators(base, 2)
        if len(item.momentum_contributions) == 2
    )

    restored = loads_completion(dumps_completion(completion))

    assert restored.momentum_contributions == completion.momentum_contributions
    assert restored.derivative_routes == completion.derivative_routes
    assert [item.derivative_degree for item in restored.momentum_contributions] == [
        2,
        2,
    ]


def test_version_one_derivative_route_is_upgraded_on_read():
    completion = completion_with_route()
    route = completion.derivative_routes[0]
    record = completion_to_record(completion)
    record["version"] = 1
    record["derivative_routes"] = [route_to_data(route)]
    del record["momentum_contributions"]

    restored = completion_from_record(record)

    assert restored.derivative_routes == (route,)
    assert restored.momentum_contributions == (
        PropagatorContribution.from_derivative_route(route),
    )


def test_completion_jsonl_round_trip(tmp_path):
    completions = [
        completion_with_route(),
        next(operator_completions(EFF_OPERATORS["2"])),
    ]
    path = tmp_path / "completions.jsonl"

    write_completion_jsonl(path, completions)
    restored = read_completion_jsonl(path)
    streamed = list(iter_completion_jsonl(path))

    assert [completion_fingerprint(item) for item in restored] == [
        completion_fingerprint(item) for item in completions
    ]
    assert [completion_fingerprint(item) for item in streamed] == [
        completion_fingerprint(item) for item in completions
    ]


def test_completion_jsonl_iterator_is_lazy(tmp_path):
    missing = tmp_path / "missing.jsonl"

    streamed = iter_completion_jsonl(missing)

    with pytest.raises(FileNotFoundError):
        next(streamed)


def test_completion_jsonl_preserves_order_and_ignores_blank_lines(tmp_path):
    first = completion_with_route()
    second = next(operator_completions(EFF_OPERATORS["2"]))
    path = tmp_path / "completions.jsonl"
    path.write_text(
        "\n  \n"
        + dumps_completion(first)
        + "\n\t\n"
        + dumps_completion(second)
        + "\n",
        encoding="utf-8",
    )

    streamed = list(iter_completion_jsonl(path))

    assert [completion_fingerprint(item) for item in streamed] == [
        completion_fingerprint(first),
        completion_fingerprint(second),
    ]


def test_completion_jsonl_reports_physical_line_for_invalid_json(tmp_path):
    path = tmp_path / "completions.jsonl"
    path.write_text("\n\n{not-json}\n", encoding="utf-8")

    with pytest.raises(
        CompletionJSONLError, match=r"completions\.jsonl:3:"
    ) as exc:
        list(iter_completion_jsonl(path))

    assert exc.value.line_number == 3
    assert isinstance(exc.value.__cause__, json.JSONDecodeError)


def test_completion_jsonl_reports_line_for_invalid_schema(tmp_path):
    path = tmp_path / "completions.jsonl"
    path.write_text(
        "\n" + json.dumps({"schema": "unknown", "version": 1}) + "\n",
        encoding="utf-8",
    )

    with pytest.raises(
        CompletionJSONLError, match=r"completions\.jsonl:2:"
    ) as exc:
        list(iter_completion_jsonl(path))

    assert "Not a neutrinomass completion record" in str(exc.value)


def test_completion_jsonl_reports_line_for_structurally_invalid_record(tmp_path):
    path = tmp_path / "completions.jsonl"
    path.write_text("\n[]\n", encoding="utf-8")

    with pytest.raises(CompletionJSONLError) as exc:
        list(iter_completion_jsonl(path))

    assert exc.value.line_number == 2
    assert isinstance(exc.value.__cause__, AttributeError)


def test_completion_jsonl_preserves_charge_validation(tmp_path):
    record = completion_to_record(completion_with_route())
    record["operator"]["tensors"][0]["charges"]["y"] = (
        "__import__('os').system('echo unsafe')"
    )
    path = tmp_path / "completions.jsonl"
    path.write_text(json.dumps(record) + "\n", encoding="utf-8")

    with pytest.raises(CompletionJSONLError, match="Invalid rational charge") as exc:
        list(iter_completion_jsonl(path))

    assert exc.value.line_number == 1


def test_completion_jsonl_streamed_read_write_preserves_fingerprints(tmp_path):
    completions = [
        completion_with_route(),
        next(operator_completions(EFF_OPERATORS["2"])),
    ]
    source = tmp_path / "source.jsonl"
    destination = tmp_path / "destination.jsonl"
    write_completion_jsonl(source, completions)

    write_completion_jsonl(destination, iter_completion_jsonl(source))

    assert [
        completion_fingerprint(item)
        for item in iter_completion_jsonl(destination)
    ] == [completion_fingerprint(item) for item in completions]
    assert destination.read_bytes() == source.read_bytes()


def test_read_completion_jsonl_remains_a_list_wrapper(tmp_path):
    path = tmp_path / "completion.jsonl"
    write_completion_jsonl(path, [completion_with_route()])

    restored = read_completion_jsonl(path)

    assert isinstance(restored, list)
    assert len(restored) == 1


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
