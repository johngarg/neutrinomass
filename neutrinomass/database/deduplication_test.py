#!/usr/bin/env python3

from copy import deepcopy

import networkx as nx
import pytest
from sympy import Rational

from neutrinomass.completions.completions import (
    are_equivalent_completions,
    expand_propagator_denominators,
    exact_completion_bucket_key,
    operator_completions,
)
from neutrinomass.completions.core import (
    Completion,
    ComplexScalar,
    DerivativeRoute,
    EffectiveOperator,
    LorentzProjection,
)
from neutrinomass.completions.fingerprints import (
    completion_fingerprint,
    completion_digest,
    lagrangian_fingerprint,
)
from neutrinomass.completions.operators import EFF_OPERATORS
from neutrinomass.database import (
    CompletionJSONLError,
    deduplicate_completion_jsonl,
    iter_completion_jsonl,
    write_completion_jsonl,
)
from neutrinomass.tensormethod import H, eps


def relabelled_scalar_completion(label, field_index, higgs_index):
    field = ComplexScalar(
        label,
        field_index,
        charges={"y": Rational(1, 2), "3b": 0},
    )
    interaction = (
        field.conj * H(higgs_index) * eps(f"-{field_index} -{higgs_index}")
    )
    return Completion(
        operator=EffectiveOperator("test", EFF_OPERATORS["1"].operator),
        partition=(),
        graph=nx.Graph(),
        exotics={field},
        terms=[interaction],
        topology="test_1",
    )


def test_disk_deduplication_preserves_first_occurrence_and_digests(tmp_path):
    first = next(operator_completions(EFF_OPERATORS["1"]))
    second = next(operator_completions(EFF_OPERATORS["2"]))
    source = tmp_path / "raw.jsonl"
    destination = tmp_path / "unique.jsonl"
    write_completion_jsonl(source, [first, first, second, first])

    report = deduplicate_completion_jsonl(source, destination, work_dir=tmp_path)
    restored = list(iter_completion_jsonl(destination))

    assert report["input_records"] == 4
    assert report["exact_classes"] == 2
    assert report["input_completion_digest"] == completion_digest(
        [first, first, second, first]
    )
    assert report["exact_completion_digest"] == completion_digest([first, second])
    assert [item.operator.name for item in restored] == [
        first.operator.name,
        second.operator.name,
    ]
    assert all(
        are_equivalent_completions(left, right)
        for left, right in zip(restored, [first, second])
    )


def test_disk_deduplication_resolves_candidate_hash_collisions_exactly(
    tmp_path, monkeypatch
):
    import neutrinomass.database.deduplication as deduplication

    left = next(operator_completions(EFF_OPERATORS["1"]))
    right = Completion(
        operator=left.operator,
        partition=left.partition,
        graph=left.graph,
        exotics=left.exotics,
        terms=[
            left.terms[0] * left.operator.fields[0].fresh_indices(),
            *left.terms[1:],
        ],
        topology=left.topology,
        canonical_topology=left.canonical_topology,
    )
    pair = (left, right)
    assert exact_completion_bucket_key(left) == exact_completion_bucket_key(right)
    assert not are_equivalent_completions(left, right)
    monkeypatch.setattr(
        deduplication,
        "lagrangian_fingerprint",
        lambda item: ("x",),
    )
    source = tmp_path / "collision.jsonl"
    destination = tmp_path / "unique.jsonl"
    write_completion_jsonl(source, pair)

    report = deduplicate_completion_jsonl(source, destination, work_dir=tmp_path)

    assert report["candidate_hash_matches"] == 1
    assert report["exact_isomorphism_comparisons"] == 1
    assert report["exact_classes"] == 2
    assert len(list(iter_completion_jsonl(destination))) == 2


def test_disk_deduplication_compares_equivalent_different_hashes(tmp_path):
    first = relabelled_scalar_completion("a", "i0", "i1")
    second = relabelled_scalar_completion("x", "i2", "i3")
    assert are_equivalent_completions(first, second)
    assert lagrangian_fingerprint(first) != lagrangian_fingerprint(second)
    source = tmp_path / "relabelled.jsonl"
    destination = tmp_path / "unique.jsonl"
    write_completion_jsonl(source, [first, second])

    report = deduplicate_completion_jsonl(source, destination, work_dir=tmp_path)

    assert report["exact_classes"] == 1
    assert report["candidate_hash_matches"] == 0
    assert report["exact_isomorphism_comparisons"] == 1


def test_disk_deduplication_does_not_collapse_different_operators(tmp_path):
    first = next(operator_completions(EFF_OPERATORS["1"]))
    renamed = deepcopy(first)
    renamed.operator = EffectiveOperator("different", first.operator.operator)
    source = tmp_path / "operators.jsonl"
    destination = tmp_path / "unique.jsonl"
    write_completion_jsonl(source, [first, renamed])

    report = deduplicate_completion_jsonl(source, destination, work_dir=tmp_path)

    assert report["exact_classes"] == 2
    assert [
        item.operator.name for item in iter_completion_jsonl(destination)
    ] == [first.operator.name, "different"]


def test_disk_deduplication_keeps_distinct_lorentz_projections(tmp_path):
    first = next(operator_completions(EFF_OPERATORS["1"]))
    second = deepcopy(first)
    first.lorentz_projection = LorentzProjection(
        ("basis-0", "basis-1"), ("1", "0"), "H", "IBP", "EOM"
    )
    second.lorentz_projection = LorentzProjection(
        ("basis-0", "basis-1"), ("0", "1"), "H", "IBP", "EOM"
    )
    source = tmp_path / "projections.jsonl"
    destination = tmp_path / "unique.jsonl"
    write_completion_jsonl(source, [first, second])

    report = deduplicate_completion_jsonl(source, destination, work_dir=tmp_path)

    assert report["exact_classes"] == 2


def test_disk_deduplication_keeps_distinct_propagator_orders(tmp_path):
    base = next(operator_completions(EFF_OPERATORS["1"]))
    first = expand_propagator_denominators(base, 1)[0]
    second = expand_propagator_denominators(base, 2)[0]
    source = tmp_path / "propagator-orders.jsonl"
    destination = tmp_path / "unique.jsonl"
    write_completion_jsonl(source, [first, second])

    assert completion_fingerprint(first) != completion_fingerprint(second)
    assert not are_equivalent_completions(first, second)

    report = deduplicate_completion_jsonl(source, destination, work_dir=tmp_path)

    assert report["exact_classes"] == 2


def test_disk_deduplication_can_collapse_local_and_routed_representatives(tmp_path):
    local = next(operator_completions(EFF_OPERATORS["1"]))
    routed = Completion(
        operator=local.operator,
        partition=local.partition,
        graph=local.graph,
        exotics=local.exotics,
        terms=local.terms,
        topology=local.topology,
        canonical_topology=local.canonical_topology,
        derivative_routes=(
            DerivativeRoute((0, 1), "x", "10", "L", "10"),
        ),
    )
    source = tmp_path / "routes.jsonl"
    destination = tmp_path / "unique.jsonl"
    write_completion_jsonl(source, [local, routed])

    report = deduplicate_completion_jsonl(source, destination, work_dir=tmp_path)

    assert report["exact_classes"] == 1


def test_disk_deduplication_cleans_temporary_files_after_malformed_input(tmp_path):
    source = tmp_path / "malformed.jsonl"
    destination = tmp_path / "unique.jsonl"
    source.write_text("{not-json}\n", encoding="utf-8")

    with pytest.raises(CompletionJSONLError):
        deduplicate_completion_jsonl(source, destination, work_dir=tmp_path)

    assert not destination.exists()
    assert sorted(path.name for path in tmp_path.iterdir()) == ["malformed.jsonl"]


def test_disk_deduplication_rejects_invalid_commit_interval(tmp_path):
    with pytest.raises(ValueError, match="commit_interval must be positive"):
        deduplicate_completion_jsonl(
            tmp_path / "source.jsonl",
            tmp_path / "destination.jsonl",
            commit_interval=0,
        )


def test_disk_deduplication_rejects_in_place_replacement(tmp_path):
    path = tmp_path / "completion.jsonl"

    with pytest.raises(ValueError, match="source and destination must differ"):
        deduplicate_completion_jsonl(path, path)
