"""Recoverable additional-mechanism searches for one admitted species multiset.

Absence is certified only after exhausting the recorded finite search domain.
An interruption, unsupported amplitude or corrupt checkpoint cannot become an
accepted candidate. Each excluded candidate has a serialized physical witness.
"""

from contextlib import closing
from functools import lru_cache
import json

from neutrinomass.completions.amplitudes import audit_amplitude_symmetrisation
from neutrinomass.completions.selection import generation_selection
from neutrinomass.completions.species import candidate_completion_stream
from neutrinomass.database.checkpoints import CheckpointStore,WorkBudget,canonical_json,identity_hash
from neutrinomass.database.global_equivalence import species_descriptors
from neutrinomass.database.planning import generation_plan
from neutrinomass.database.rebuild import operator_registry,operator_scale_tev,validate_completion
from neutrinomass.database.serialization import dumps_completion,loads_completion,tensor_to_data


def inventory(candidate):
    return sorted((tensor_to_data(f) for f in candidate.exotics),key=canonical_json)


@lru_cache(maxsize=256)
def _plan(name):
    ranges=generation_plan(name)
    return ranges,dict(units=len(ranges),sha256=identity_hash(ranges))


def checked_tree_witness(completion,available,operator):
    decoded=loads_completion(dumps_completion(completion))
    validate_completion(decoded,check_vanishing=False)
    if any(term.safe_simplify()==0 for term in decoded.terms):
        return None
    audit=audit_amplitude_symmetrisation(decoded)
    if audit.status=="unsupported":
        raise ValueError(f"UNSUPPORTED_CANDIDATE_MECHANISM: {audit.reason}")
    if audit.is_zero:
        return None
    required=species_descriptors(decoded)
    if any(available.get(label)!=descriptor for label,descriptor in required.items()):
        raise ValueError("CANDIDATE_MECHANISM_CHANGED_SPECIES")
    return dict(code="HIGHER_SCALE_CANDIDATE_TREE",verified=True,operator=operator.name,
        scale_tev=operator_scale_tev(operator),
        required_species=tuple(sorted(required.values(),key=repr)),
        species_assignment=dict(sorted(required.items())),
        record=dumps_completion(decoded),amplitude_witness=audit.witness)


def tree_search(candidate,higher_operators,directory,contract,*,budget=None,plans=None):
    """Continue at a completed assignment, including assignments with no output.

    ``plans`` can supply explicit bounded controls; their exact ranges are part
    of the identity and cannot be reused as an exhaustive search of the domain.
    Production omits it and uses every ordinary generation range.
    """
    budget=budget or WorkBudget()
    higher_operators=tuple(higher_operators)
    operators=operator_registry()
    selections={name:_plan(name)[0] for name in higher_operators} if plans is None else plans
    if set(selections)!=set(higher_operators):
        raise ValueError("MECHANISM_PLAN_INVENTORY_MISMATCH")
    identity=dict(contract,algorithm="candidate-tree-search-v1",inventory=inventory(candidate),
        higher_operators=higher_operators,
        plans={name:(_plan(name)[1] if plans is None else
                     dict(units=len(selections[name]),sha256=identity_hash(selections[name])))
               for name in higher_operators},exhaustive=plans is None)
    available=species_descriptors(candidate)
    with CheckpointStore(directory,identity) as store:
        if store.state["complete"]:
            return store.state
        progress=dict(store.state["progress"] or dict(operator_index=0,unit_index=0,
            next_assignment=None,generated=0,witness=False))
        def save():
            store.save(progress)
        for operator_index in range(progress["operator_index"],len(higher_operators)):
            name=higher_operators[operator_index]
            for unit_index in range(progress["unit_index"],len(selections[name])):
                selection=dict(selections[name][unit_index])
                if progress["next_assignment"] is not None:
                    selection["assignment_start"]=progress["next_assignment"]
                def boundary(partition):
                    progress["next_assignment"]=partition["assignment"]+1
                    budget.boundary(save,progress=dict(stage="candidate-tree",operator=name,
                        **progress))
                with generation_selection(selection,boundary),closing(
                        candidate_completion_stream(operators[name],candidate.exotics)) as stream:
                    for completion in stream:
                        progress["generated"]+=1
                        witness=checked_tree_witness(completion,available,operators[name])
                        if witness is not None:
                            artifact=store.publish_lines([canonical_json(witness)])
                            progress["witness"]=True
                            store.save(progress,artifacts={"witness":artifact},complete=True)
                            return store.state
                progress.update(unit_index=unit_index+1,next_assignment=None)
                budget.boundary(save,progress=dict(stage="candidate-tree",operator=name,**progress))
            progress.update(operator_index=operator_index+1,unit_index=0,next_assignment=None)
        store.save(progress,complete=True)
        return store.state


def loop_search(candidate,directory,contract,*,budget=None):
    from neutrinomass.database.loop_mechanisms import candidate_loop_witness
    budget=budget or WorkBudget()
    identity=dict(contract,algorithm="candidate-heavy-loop-search-v2",inventory=inventory(candidate))
    with CheckpointStore(directory,identity) as store:
        if store.state["complete"]:
            return store.state
        progress=dict(store.state["progress"] or dict(next_cycle=0,witness=False))
        def save():
            store.save(progress)
        def boundary(next_cycle):
            progress["next_cycle"]=next_cycle
            budget.boundary(save,progress=dict(stage="candidate-heavy-loop",**progress))
        witness=candidate_loop_witness(candidate,start_cycle=progress["next_cycle"],boundary=boundary)
        artifacts={}
        if witness is not None:
            artifacts["witness"]=store.publish_lines([canonical_json(witness)])
            progress["witness"]=True
        store.save(progress,artifacts=artifacts,complete=True)
        return store.state


def read_witness(directory,state):
    from pathlib import Path
    from neutrinomass.database.resumable_reduction import lines_from_paths
    if not state["complete"]:
        raise ValueError("MECHANISM_SEARCH_INCOMPLETE")
    if not state["progress"]["witness"]:
        return None
    rows=list(lines_from_paths([Path(directory)/state["artifacts"]["witness"]["path"]]))
    if len(rows)!=1:
        raise ValueError("MECHANISM_WITNESS_COUNT_MISMATCH")
    return json.loads(rows[0])
