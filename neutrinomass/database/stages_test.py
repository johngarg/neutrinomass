import pytest

from neutrinomass.completions.completions import get_topology_data
from neutrinomass.completions.operators import EFF_OPERATORS
from neutrinomass.database.checkpoints import BudgetExpired, WorkBudget
from neutrinomass.database.stages import generate_unit, map_records, output_paths, decode_record, amplitude_record
from neutrinomass.database.resumable_reduction import lines_from_paths


def test_generation_interrupt_and_decode_match_uninterrupted(tmp_path):
    topology = get_topology_data(**EFF_OPERATORS["1"].topology_type)[0]["canonical_topology"]
    selection = dict(topology=topology,colour=0,assignment_start=0,assignment_stop=None)
    full = generate_unit("1",selection,tmp_path/"full",{})
    with pytest.raises(BudgetExpired):
        generate_unit("1",selection,tmp_path/"resumed",{},budget=WorkBudget(fault_after=1))
    resumed = generate_unit("1",selection,tmp_path/"resumed",{})
    full_lines = list(lines_from_paths(output_paths(tmp_path/"full",full)))
    assert full_lines == list(lines_from_paths(output_paths(tmp_path/"resumed",resumed)))
    decoded = map_records(output_paths(tmp_path/"resumed",resumed),tmp_path/"decode",{},
                          decode_record,algorithm="decoded-validation-v1")
    assert decoded["progress"]["processed"] == len(full_lines)
    assert decoded["progress"]["accepted"] > 0


def test_assignment_ranges_reconstruct_ordered_generation(tmp_path):
    topology = get_topology_data(**EFF_OPERATORS["1"].topology_type)[0]["canonical_topology"]
    base = dict(topology=topology,colour=0,assignment_start=0,assignment_stop=None)
    full = generate_unit("1",base,tmp_path/"full",{})
    pieces=[]
    for i,(start,stop) in enumerate([(0,1),(1,None)]):
        state = generate_unit("1",dict(base,assignment_start=start,assignment_stop=stop),tmp_path/str(i),{})
        pieces.extend(lines_from_paths(output_paths(tmp_path/str(i),state)))
    assert pieces == list(lines_from_paths(output_paths(tmp_path/"full",full)))
