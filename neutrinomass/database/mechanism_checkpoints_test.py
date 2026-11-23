from pathlib import Path
from types import SimpleNamespace
import json

import pytest

from neutrinomass.completions.core import ComplexScalar
from neutrinomass.database.checkpoints import BudgetExpired,WorkBudget
from neutrinomass.database.mechanism_checkpoints import tree_search,loop_search,read_witness
from neutrinomass.database.planning import generation_plan
from neutrinomass.database.serialization import iter_completion_jsonl,loads_completion
from neutrinomass.completions.amplitudes import audit_amplitude_symmetrisation
from neutrinomass.database.filtering import multiset_contains
from neutrinomass.database.global_equivalence import species_descriptors
from neutrinomass.database.loop_mechanisms_test import scotogenic


def test_absence_search_resumes_assignments_that_produce_no_records(tmp_path):
    candidate=SimpleNamespace(exotics={ComplexScalar("absent","",charges={"y":7,"3b":0})})
    whole=tree_search(candidate,["1"],tmp_path/"whole",{})
    assert whole["complete"] and not whole["progress"]["witness"]
    with pytest.raises(BudgetExpired):
        tree_search(candidate,["1"],tmp_path/"resume",{},budget=WorkBudget(fault_after=2))
    interrupted=json.loads((tmp_path/"resume/checkpoint.json").read_text())
    assert not interrupted["complete"]
    assert interrupted["progress"]["next_assignment"]>0
    resumed=tree_search(candidate,["1"],tmp_path/"resume",{})
    assert resumed["progress"]==whole["progress"]
    assert read_witness(tmp_path/"resume",resumed) is None
    with pytest.raises(ValueError,match="Incompatible"):
        tree_search(candidate,["3a"],tmp_path/"resume",{})


def test_copy_witness_survives_persistent_search_and_serialization(tmp_path):
    fixture=Path(__file__).resolve().parents[2]/"tests/fixtures/copy-mechanisms/0.jsonl"
    candidate=next(iter_completion_jsonl(fixture));name=candidate.operator.name
    plans={name:[s for s in generation_plan(name) if s["topology"]==candidate.canonical_topology]}
    state=tree_search(candidate,[name],tmp_path,{},plans=plans)
    witness=read_witness(tmp_path,state)
    decoded=loads_completion(witness["record"])
    assert audit_amplitude_symmetrisation(decoded).status=="nonzero"
    assert state["contract"]["exhaustive"] is False
    assert multiset_contains(tuple(species_descriptors(candidate).values()),
                             tuple(species_descriptors(decoded).values()))


def test_loop_witness_publishes_once_and_rejects_incompatible_candidate(tmp_path):
    result=loop_search(scotogenic(),tmp_path,{})
    witness=read_witness(tmp_path,result)
    assert witness["verified"] and witness["colour_dimension"]==1
    assert loop_search(scotogenic(),tmp_path,{})==result
    with pytest.raises(ValueError,match="Incompatible"):
        loop_search(scotogenic(baryon=3),tmp_path,{})
