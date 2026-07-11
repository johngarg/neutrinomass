import pytest

from neutrinomass.completions.completions import operator_completions
from neutrinomass.completions.core import ComplexScalar, VectorLikeDiracFermion
from neutrinomass.completions.operators import EFF_OPERATORS
from neutrinomass.tensormethod import L, eps
from neutrinomass.database.deduplication import deduplicate_completion_jsonl
from neutrinomass.database.rebuild import (
    audit_exact_artifact,
    classify_historical_classes,
    file_sha256,
    filter_democratic_registry,
    operator_inventory,
    operator_registry,
    partition_bucketable_historical_records,
    write_generated_artifact,
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


def test_streamed_generation_and_disk_backed_historical_audit(tmp_path):
    raw_path = tmp_path / "raw.jsonl"
    exact_path = tmp_path / "exact.jsonl"
    model_path = tmp_path / "models.jsonl"
    completions = list(operator_completions(EFF_OPERATORS["1"]))

    generated = write_generated_artifact(raw_path, iter(completions))
    deduplication = deduplicate_completion_jsonl(
        raw_path, exact_path, work_dir=tmp_path
    )
    audit = audit_exact_artifact(
        exact_path,
        model_path,
        [completions[0]],
        work_dir=tmp_path,
    )

    assert generated["records"] == 8
    assert generated["local"] == 8
    assert generated["routed"] == 0
    assert deduplication["exact_classes"] == 3
    assert audit["records"] == 3
    assert audit["democratic_models"] == 3
    assert audit["missing_historical_fingerprints"] == []
    assert model_path.read_text(encoding="utf-8").count("\n") == 3


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

    valid, invalid = classify_historical_classes([completion])

    assert valid == []
    assert len(invalid) == 1
    assert invalid[0]["reason"] == "vanishing UV interaction"
    assert invalid[0]["vanishing_term_indices"] == [0]


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
