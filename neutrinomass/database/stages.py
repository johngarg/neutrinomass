"""Persistent generation, validation and amplitude work units.

Generation and decoding intentionally remain separate passes: constructing a
decoded SymPy tensor can affect its global tensor-head cache. A generation
unit has one producer/placement/topology/colour and a half-open assignment
range. Exact reduction across all units remains mandatory.
"""

import json
from pathlib import Path
import tempfile

from neutrinomass.database.checkpoints import CheckpointStore, WorkBudget, canonical_json
from neutrinomass.database.resumable_reduction import input_contract, lines_from_paths
from neutrinomass.database.serialization import dumps_completion, loads_completion
from neutrinomass.database.rebuild import round_trip_signature


def output_paths(directory, state, prefix="chunk-"):
    return [Path(directory)/a["path"] for k,a in sorted(state["artifacts"].items())
            if k.startswith(prefix)]


def generate_unit(operator_name, selection, directory, contract, *, budget=None):
    from neutrinomass.completions.selection import generation_selection
    from neutrinomass.completions.operators import DERIV_EFF_OPERATORS
    from neutrinomass.database.rebuild import completion_stream, validate_completion
    if selection.get("topology") is None or selection.get("colour") is None:
        raise ValueError("Generation units require a fixed topology and colour structure")
    if operator_name in DERIV_EFF_OPERATORS and (
        selection.get("producer") not in {"local","routed"}
        or (selection["producer"] == "local" and selection.get("placement") is None)
    ):
        raise ValueError("Derivative units require a fixed producer and local placement")
    budget = budget or WorkBudget()
    identity = dict(contract, algorithm="generation-v1", operator=operator_name, selection=selection)
    with CheckpointStore(directory,identity) as store:
        if store.state["complete"]:
            return store.state
        progress = dict(store.state["progress"] or dict(
            next_assignment=selection.get("assignment_start",0), records=0, chunks=0))
        resumed_selection = dict(selection, assignment_start=progress["next_assignment"])
        with tempfile.TemporaryFile(mode="w+t") as chunk:
            buffered = 0

            def save(complete=False):
                nonlocal buffered
                artifacts = {}
                if buffered:
                    chunk.seek(0)
                    artifacts[f"chunk-{progress['chunks']:08d}"] = store.publish_lines(chunk)
                    progress["chunks"] += 1
                    chunk.seek(0); chunk.truncate(); buffered = 0
                store.save(progress, artifacts=artifacts, complete=complete)

            def boundary(partition):
                progress["next_assignment"] = partition["assignment"]+1
                budget.boundary(save, progress=dict(stage="generation", **progress))

            with generation_selection(resumed_selection,boundary):
                for completion in completion_stream(operator_name):
                    validate_completion(completion, check_vanishing=False)
                    chunk.write(canonical_json(dict(
                        payload=dumps_completion(completion),
                        signature=repr(round_trip_signature(completion))))+"\n")
                    buffered += 1
                    progress["records"] += 1
            save(complete=True)
            return store.state


def map_records(paths, directory, contract, transform, *, algorithm, budget=None, scratch=None):
    """Checkpoint an ordered per-record transform and its scientific witnesses."""
    paths = tuple(Path(p) for p in paths)
    budget = budget or WorkBudget()
    identity = dict(contract, algorithm=algorithm, inputs=input_contract(paths))
    with CheckpointStore(directory,identity) as store:
        if store.state["complete"]:
            return store.state
        with store.database(scratch) as db:
            db.executescript("""
                CREATE TABLE IF NOT EXISTS output (input INTEGER, part INTEGER, payload TEXT,
                    PRIMARY KEY(input,part));
                CREATE TABLE IF NOT EXISTS evidence (input INTEGER PRIMARY KEY, payload TEXT);
            """)
            progress = dict(store.state["progress"] or dict(processed=0, accepted=0, rejected=0))
            def save():
                store.snapshot(db,progress)
            for ordinal,payload in enumerate(lines_from_paths(paths),1):
                if ordinal <= progress["processed"]:
                    continue
                output, evidence = transform(payload)
                for part,line in enumerate(output):
                    db.execute("INSERT INTO output VALUES (?,?,?)",(ordinal,part,line))
                db.execute("INSERT INTO evidence VALUES (?,?)",(ordinal,canonical_json(evidence)))
                progress["processed"] = ordinal
                progress["accepted"] += len(output)
                progress["rejected"] += not bool(output)
                budget.boundary(save, progress=dict(stage=algorithm, **progress))
            artifacts = {name:store.publish_lines(row[0] for row in db.execute(query))
                for name,query in {
                    "output":"SELECT payload FROM output ORDER BY input,part",
                    "evidence":"SELECT payload FROM evidence ORDER BY input",
                }.items()}
            store.save(progress,artifacts=artifacts)
            store.snapshot(db,progress,complete=True)
            return store.state


def decode_record(line):
    from neutrinomass.database.rebuild import validate_completion, topology_key
    record = json.loads(line)
    completion = loads_completion(record["payload"])
    signature=round_trip_signature(completion)
    if repr(signature) != record["signature"]:
        raise ValueError("SERIALIZATION_SIGNATURE_CHANGED")
    # Structural/routing errors block; the known tensor-statistics correction
    # is retained as a concrete rejected-record witness.
    validate_completion(completion,check_vanishing=False)
    zero = [i for i,t in enumerate(completion.terms) if t.safe_simplify() == 0]
    evidence = dict(code="DECODED_VANISHING_VERTEX" if zero else "ACCEPTED",
                    topology=topology_key(completion), vanishing_term_indices=zero,
                    fingerprint=repr(signature[0]),signature=repr(signature),
                    routed=bool(completion.momentum_contributions))
    if zero:
        evidence["record"] = record["payload"]
        return [],evidence
    return [dumps_completion(completion)],evidence


def amplitude_record(line):
    from neutrinomass.completions.amplitudes import audit_amplitude_symmetrisation
    completion = loads_completion(line)
    audit = audit_amplitude_symmetrisation(completion)
    if audit.status == "unsupported":
        raise ValueError(f"UNSUPPORTED_AMPLITUDE: {audit.reason}")
    evidence = dict(code="AMPLITUDE_ZERO" if audit.is_zero else "AMPLITUDE_NONZERO",
                    derivative_sectors=audit.derivative_sectors, witness=audit.witness)
    if audit.is_zero:
        evidence["record"] = line
    return ([] if audit.is_zero else [line]),evidence


def historical_transform(operator_name):
    """Only for hash-verified package-controlled historical archives."""
    import ast
    from neutrinomass.database.database import LazyCompletion
    from neutrinomass.database.rebuild import historical_record_normalizer, validate_completion, _invalid_historical_class, HistoricalRecordRejection
    from neutrinomass.completions.completions import exact_completion_bucket_key
    from neutrinomass.completions.routing import InvalidMomentumGraph
    normalize=historical_record_normalizer(operator_name)
    def transform(line):
        call=ast.parse(line,mode="eval").body
        if not isinstance(call,ast.Call) or not isinstance(call.func,ast.Name) or call.func.id != "LazyCompletion" or call.args:
            raise ValueError("Unrecognised legacy archive record")
        values={kw.arg:ast.literal_eval(kw.value) for kw in call.keywords}
        if set(values) != {"head","tail"} or values["head"]["operator_name"] != operator_name:
            raise ValueError("Legacy archive header differs from selected operator")
        item=LazyCompletion(**values)
        records=[]
        invalid=[]
        try:
            for completion in normalize(item.force(trusted=True)):
                try:
                    exact_completion_bucket_key(completion)
                except ValueError:
                    try:
                        validate_completion(completion)
                    except InvalidMomentumGraph:
                        raise
                    except ValueError as error:
                        invalid.append(_invalid_historical_class(completion,error))
                    else:
                        raise
                else:
                    records.append(dumps_completion(completion))
        except HistoricalRecordRejection as error:
            invalid.append(_invalid_historical_class(error.completion,error))
        return records,dict(code="HISTORICAL_NORMALIZED",normalized_records=len(records),
                            invalid=invalid,representation_set=sorted(set(item.head["quantum_numbers"])))
    return transform
