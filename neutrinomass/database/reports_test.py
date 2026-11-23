import pytest
from neutrinomass.database.reports import normalize_census_report, stable_diagnostic_record


def test_old_mislabelled_tev_is_migrated_without_factor_of_thousand():
    old = dict(schema_version=1,operator_scale_gev=0.8,historical={})
    new = normalize_census_report(old)
    assert new["operator_scale_tev"] == 0.8
    assert new["scale_unit"] == "TeV" and new["schema_version"] == 2
    assert old["schema_version"] == 1
    with pytest.raises(ValueError):
        normalize_census_report(dict(old,operator_scale_tev=800))


def test_only_known_label_dependent_diagnostic_is_normalised():
    a = dict(fingerprint="fixed",reason="Inconsistent quantum numbers for exotic species a")
    b = dict(a,reason="Inconsistent quantum numbers for exotic species b")
    assert stable_diagnostic_record(a) == stable_diagnostic_record(b)
    assert stable_diagnostic_record(dict(a,reason="unknown first")) != stable_diagnostic_record(dict(b,reason="unknown second"))
