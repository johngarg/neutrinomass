from copy import deepcopy
from types import SimpleNamespace

import pytest

from neutrinomass.completions.completions import operator_completions
from neutrinomass.completions.core import (
    ComplexScalar,
    VectorLikeDiracFermion,
    cons_completion_field,
)
from neutrinomass.completions.operators import DERIV_EFF_OPERATORS, EFF_OPERATORS
from neutrinomass.tensormethod import L, eps
from neutrinomass.tensormethod.core import FERMI, IndexedField
from neutrinomass.database.deduplication import deduplicate_completion_jsonl
from neutrinomass.database.export import export_completion
from neutrinomass.database.serialization import (
    dumps_completion,
    iter_completion_jsonl,
    loads_completion,
)
from neutrinomass.database.rebuild import (
    audit_exact_artifact,
    classify_historical_classes,
    completion_stream,
    file_sha256,
    filter_democratic_registry,
    operator_inventory,
    operator_registry,
    partition_bucketable_historical_records,
    write_generated_artifact,
    write_historical_artifact,
)


def test_operator_registry_covers_the_complete_packaged_catalogue():
    registry = operator_registry()

    assert len(registry) == 243
    assert len({name for name in registry if name.startswith("D")}) == 51
    assert len({name for name in registry if not name.startswith("D")}) == 192


def test_operator_inventory_requires_an_exact_legacy_mapping(tmp_path, monkeypatch):
    operator = EFF_OPERATORS["1"]
    legacy = tmp_path / "op_1.dat"
    legacy.write_bytes(b"trusted legacy placeholder")
    monkeypatch.setattr(
        "neutrinomass.database.rebuild.EFF_OPERATORS", {"1": operator}
    )
    monkeypatch.setattr(
        "neutrinomass.database.rebuild.DERIV_EFF_OPERATORS", {}
    )

    inventory = operator_inventory(tmp_path)

    assert inventory == (
        {
            "operator": "1",
            "kind": "regular",
            "derivatives": 0,
            "legacy_path": str(legacy.resolve()),
            "legacy_sha256": file_sha256(legacy),
        },
    )
    (tmp_path / "op_missing.dat").write_bytes(b"missing definition")
    with pytest.raises(ValueError) as error:
        operator_inventory(tmp_path)
    assert "missing definitions=['missing']" in str(error.value)


def test_derivative_census_uses_safe_canonical_partition_preflight(monkeypatch):
    calls = []

    def generate(operator, *, canonical_partitions):
        calls.append((operator.name, canonical_partitions))
        return ()

    monkeypatch.setattr(
        "neutrinomass.database.rebuild.deriv_operator_completion_stream", generate
    )
    monkeypatch.setattr(
        "neutrinomass.database.rebuild.DERIV_EFF_OPERATORS",
        {
            "verified": DERIV_EFF_OPERATORS["D20"],
            "one": DERIV_EFF_OPERATORS["D5b"],
            "projected": DERIV_EFF_OPERATORS["D8g"],
            "several": DERIV_EFF_OPERATORS["D21"],
        },
    )

    assert list(completion_stream("verified")) == []
    assert list(completion_stream("one")) == []
    assert list(completion_stream("projected")) == []
    assert list(completion_stream("several")) == []
    assert calls == [
        ("D20", True),
        ("D5b", False),
        ("D8g", False),
        ("D21", False),
    ]


def test_streamed_generation_and_disk_backed_historical_audit(
    tmp_path, monkeypatch
):
    raw_path = tmp_path / "raw.jsonl"
    exact_path = tmp_path / "exact.jsonl"
    model_path = tmp_path / "models.jsonl"
    completions = list(operator_completions(EFF_OPERATORS["1"]))

    generated = write_generated_artifact(raw_path, iter(completions))
    deduplication = deduplicate_completion_jsonl(
        raw_path, exact_path, work_dir=tmp_path
    )
    monkeypatch.setattr(
        "neutrinomass.database.rebuild.audit_amplitude_symmetrisation",
        lambda completion: pytest.fail(
            "matched historical classes inherit the audited exact survivor"
        ),
    )
    audit = audit_exact_artifact(
        exact_path,
        model_path,
        [completions[0]],
        work_dir=tmp_path,
        classify_historical=True,
    )

    assert generated["records"] == 8
    assert generated["local"] == 8
    assert generated["routed"] == 0
    assert deduplication["exact_classes"] == 3
    assert audit["records"] == 3
    assert audit["democratic_models"] == 3
    assert audit["valid_historical_classes"] == 1
    assert audit["reproduced_historical_classes"] == 1
    assert audit["missing_historical_fingerprints"] == []
    assert model_path.read_text(encoding="utf-8").count("\n") == 3


def test_legacy_historical_normalisation_is_streamed_to_safe_jsonl(tmp_path):
    completion = next(operator_completions(EFF_OPERATORS["1"]))
    legacy = tmp_path / "op_1.dat"
    legacy.write_text(
        (export_completion(completion) + "\n") * 2,
        encoding="utf-8",
    )
    destination = tmp_path / "historical.jsonl"

    report = write_historical_artifact("1", legacy, destination)

    assert report["records"] == 2
    assert report["bucketable_records"] == 2
    assert report["unbucketable_invalid"] == []
    assert len(list(iter_completion_jsonl(destination))) == 2


def test_generated_artifact_rejects_vertices_zero_after_round_trip(
    tmp_path, monkeypatch
):
    valid = next(operator_completions(EFF_OPERATORS["1"]))
    invalid = deepcopy(valid)
    scalar = ComplexScalar("phi", "-c0 i0", charges={"y": 0, "3b": 0})
    fermion = VectorLikeDiracFermion(
        "psi", "u1 -c2 -c1", charges={"y": 0, "3b": 0}
    )
    invalid.terms = [
        scalar
        * L("u0 i1")
        * fermion
        * eps("-u0 -u1")
        * eps("-i0 -i1")
        * eps("c0 c1 c2")
    ]
    path = tmp_path / "generated.jsonl"
    monkeypatch.setattr(
        "neutrinomass.database.rebuild.validate_completion",
        lambda completion, check_vanishing=True: None,
    )

    report = write_generated_artifact(
        path, [valid, invalid], terms_prevalidated=True
    )

    assert report["source_records"] == 2
    assert report["records"] == 1
    assert report["decoded_vertex_rejections"] == 1
    assert report["decoded_rejection_topologies"]
    survivor = list(iter_completion_jsonl(path))[0]
    assert all(term.safe_simplify() != 0 for term in survivor.terms)


def test_safe_round_trip_corrects_legacy_exotic_tensor_symmetry():
    completion = next(operator_completions(EFF_OPERATORS["1"]))
    charges = {"y": 1, "3b": 0}
    legacy_kwargs = {
        "symmetry": [[1], [1], [1, 1]],
        "charges": charges,
        "nf": 1,
        "dynkin": "01012",
        "comm": FERMI,
        "latex": "psi",
        "is_conj": True,
        "is_unbarred": True,
    }
    first = cons_completion_field(
        IndexedField(
            label="legacy_psi†",
            indices="d0 -c0 i0 i1",
            **legacy_kwargs,
        )
    )
    second = cons_completion_field(
        IndexedField(
            label="legacy_psi†",
            indices="d1 -c1 i2 i3",
            **legacy_kwargs,
        )
    )
    eta = ComplexScalar(
        "legacy_eta", "-c2", charges={"y": -2, "3b": 0}
    )
    legacy_term = (
        eta
        * first
        * second
        * eps("-d0 -d1")
        * eps("-i2 -i1")
        * eps("-i3 -i0")
        * eps("c0 c1 c2")
    )
    completion.exotics = {eta, first, second}
    completion.terms = [legacy_term]

    restored = loads_completion(dumps_completion(completion))

    assert legacy_term.safe_simplify() != 0
    assert restored.terms[0].safe_simplify() == 0


def test_regenerated_democratic_filter_uses_surviving_upstream_subsets():
    registry = {
        "high": {("F",): (("top", "top"),)},
        "equal": {("F", "E"): (("top", "top"),)},
        "low": {
            ("F", "S"): (("top", "top"),),
            ("W",): (("top", "top"),),
            ("X",): (("top", "top"),),
        },
    }
    scales = {"high": 10.0, "equal": 10.0, "low": 1.0}

    filtered = filter_democratic_registry(registry, scales, [("W",)])

    assert set(filtered["survivors"]["high"]) == {("F",)}
    assert set(filtered["survivors"]["equal"]) == {("F", "E")}
    assert set(filtered["survivors"]["low"]) == {("X",)}
    assert filtered["removed_by_mass"] == {"low": 1}
    assert filtered["removed_by_one_loop_weinberg"] == {"low": 1}


def test_historical_audit_classifies_vanishing_legacy_classes():
    scalar = ComplexScalar("phi", "-c0 i0", charges={"y": 0, "3b": 0})
    fermion = VectorLikeDiracFermion(
        "psi", "u1 -c2 -c1", charges={"y": 0, "3b": 0}
    )
    vanishing = (
        scalar
        * L("u0 i1")
        * fermion
        * eps("-u0 -u1")
        * eps("-i0 -i1")
        * eps("c0 c1 c2")
    )
    completion = next(operator_completions(EFF_OPERATORS["1"]))
    completion.terms = [vanishing]

    valid, invalid, unsupported = classify_historical_classes([completion])

    assert valid == []
    assert len(invalid) == 1
    assert invalid[0]["reason"] == "vanishing UV interaction"
    assert invalid[0]["vanishing_term_indices"] == [0]
    assert unsupported == []


def test_historical_audit_defers_unsupported_provenance_to_coverage(monkeypatch):
    completion = next(operator_completions(EFF_OPERATORS["1"]))
    monkeypatch.setattr(
        "neutrinomass.database.rebuild.audit_amplitude_symmetrisation",
        lambda completion: SimpleNamespace(
            status="unsupported", is_zero=False, reason="legacy provenance"
        ),
    )

    valid, invalid, unsupported = classify_historical_classes([completion])

    assert valid == [completion]
    assert invalid == []
    assert len(unsupported) == 1
    assert unsupported[0]["reason"] == "legacy provenance"


def test_historical_audit_classifies_inconsistent_exotic_labels():
    completions = list(operator_completions(EFF_OPERATORS["1"]))
    completion = completions[4]
    completion.exotics = {
        next(iter(completions[4].exotics)),
        next(iter(completions[5].exotics)),
    }
    completion.terms = [completions[4].terms[0], completions[5].terms[0]]

    bucketable, invalid = partition_bucketable_historical_records([completion])

    assert bucketable == []
    assert len(invalid) == 1
    assert invalid[0]["reason"] == (
        "Inconsistent quantum numbers for exotic species ψ"
    )
    assert invalid[0]["vanishing_term_indices"] == []
