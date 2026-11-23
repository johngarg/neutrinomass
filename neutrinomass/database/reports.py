"""Explicit migration of census report units and deterministic diagnostics."""

from copy import deepcopy


def diagnostic(reason):
    reason = str(reason)
    if reason.startswith("Inconsistent quantum numbers for exotic species "):
        return {"code":"inconsistent_species_quantum_numbers"}
    known = {
        "vanishing UV interaction":"vanishing_uv_interaction",
        "non-renormalisable or non-singlet UV interaction":"invalid_uv_vertex",
    }
    if reason in known:
        return {"code":known[reason]}
    # Unknown errors remain fully significant. Never discard arbitrary text.
    return {"code":"unclassified_validation_error","detail":reason}


def normalize_census_report(report):
    report = deepcopy(report)
    version = report.get("schema_version")
    if version == 1:
        if "operator_scale_tev" in report:
            raise ValueError("Ambiguous units in version-one report")
        # Historical census schema 1 labelled TeV numerical values as GeV.
        report["operator_scale_tev"] = report.pop("operator_scale_gev")
        report["schema_version"] = 2
        report["scale_unit"] = "TeV"
        report["migration"] = {"source_schema":1,"operation":"correct_tev_label_without_rescaling"}
    elif version != 2:
        raise ValueError("Unsupported census report schema")
    if report.get("scale_unit") != "TeV" or "operator_scale_gev" in report:
        raise ValueError("Census report does not have explicit unambiguous TeV units")
    for record in report.get("historical",{}).get("invalid",[]):
        record.setdefault("diagnostic",diagnostic(record["reason"]))
    return report


def stable_diagnostic_record(record):
    record = deepcopy(record)
    record.setdefault("diagnostic",diagnostic(record["reason"]))
    record.pop("reason")
    return record
