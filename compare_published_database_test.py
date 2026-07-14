from compare_published_database import build_operator_rows, build_report


def test_comparison_reports_common_operator_deltas_and_current_only_rows():
    census = {
        "historical": {"records": 5},
        "structural_exact_classes": 8,
        "amplitude_symmetrisation": {"rejected_records": 2},
        "exact_classes": {"all": 6},
        "species_models": 5,
        "democratic_models": 4,
        "generation_rejections": {"decoded_vanishing_uv_interactions": 1},
    }
    manifest = {
        "source_commit": "science",
        "operators": {
            "D3": {"kind": "derivative", "derivatives": 1, "census": census},
            "D13c": {
                "kind": "derivative",
                "derivatives": 1,
                "census": census,
            },
        },
    }
    filtering = {
        "totals": {"input": 8, "survivors": 3},
        "operators": {"D3": {"survivors": 2}, "D13c": {"survivors": 1}},
    }

    rows = build_operator_rows(manifest, filtering, {"D3": 3})
    report = build_report(
        manifest,
        filtering,
        rows,
        {
            "physical_lagrangians": 430811,
            "physical_lagrangian_operator_incidences": 9,
            "democratic_models": 141990,
            "filtered_lagrangians": 11484,
            "filtered_lagrangian_operator_incidences": 3,
            "filtered_democratic_models": 11217,
        },
    )

    by_name = {row["operator"]: row for row in rows}
    assert by_name["D3"]["model_delta"] == -1
    assert by_name["D3"]["filtered_delta"] == -1
    assert by_name["D13c"]["published_models"] is None
    assert by_name["D13c"]["filtered_delta"] is None
    assert report["headline_delta"] == {
        "physical_lagrangians": 1,
        "democratic_models": 1,
        "filtered_lagrangians": 1,
        "filtered_democratic_models": 1,
    }
