from copy import deepcopy
import pytest
from scripts.validate_71p_gate import report_gate,GateFailure


def report():
    return dict(operator="71p",hash_seed="0",source_commit="tested",records=dict(generator=758208),
        stage_manifest=dict(generation_units=3168,topologies=[str(i) for i in range(264)]),
        historical=dict(missing=0,amplitude_unsupported=[],valid_classes=2,reproduced=2),
        exact_classes=dict(all=3))


@pytest.mark.parametrize("mutation",[
    lambda r:r["stage_manifest"]["topologies"].pop(),
    lambda r:r["records"].update(generator=758207),
    lambda r:r["historical"].update(missing=1),
    lambda r:r["historical"].update(amplitude_unsupported=["unresolved"]),
    lambda r:r.update(hash_seed="1"),
    lambda r:r.update(source_commit="untested"),
])
def test_incomplete_or_unexplained_pilot_cannot_open_production(mutation):
    candidate=report()
    report_gate(candidate,seed=0,producer="tested",unit_count=3168)
    mutation(candidate)
    with pytest.raises(GateFailure):
        report_gate(candidate,seed=0,producer="tested",unit_count=3168)
