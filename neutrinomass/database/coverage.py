"""Historical coverage with a durable cursor inside exact-match buckets."""

from pathlib import Path

from neutrinomass.completions.completions import are_equivalent_completions, exact_completion_bucket_key
from neutrinomass.completions.amplitudes import audit_amplitude_symmetrisation
from neutrinomass.completions.routing import InvalidMomentumGraph
from neutrinomass.database.checkpoints import CheckpointStore, WorkBudget, canonical_json
from neutrinomass.database.resumable_reduction import lines_from_paths, input_contract
from neutrinomass.database.global_equivalence import DecodedCache
from neutrinomass.database.serialization import loads_completion
from neutrinomass.database.rebuild import validate_completion, _invalid_historical_class


def historical_coverage(exact_paths, historical_paths, directory, contract, *, budget=None, scratch=None):
    exact_paths, historical_paths = tuple(exact_paths), tuple(historical_paths)
    budget = budget or WorkBudget()
    identity = dict(contract,algorithm="historical-coverage-v1",
                    exact=input_contract(exact_paths),historical=input_contract(historical_paths))
    with CheckpointStore(directory,identity) as store:
        if store.state["complete"]:
            return store.state
        with store.database(scratch) as db:
            db.executescript("""
                CREATE TABLE IF NOT EXISTS exact (id INTEGER PRIMARY KEY,bucket TEXT,payload TEXT);
                CREATE INDEX IF NOT EXISTS exact_buckets ON exact(bucket,id);
                CREATE TABLE IF NOT EXISTS coverage (id INTEGER PRIMARY KEY,status TEXT,evidence TEXT);
            """)
            progress=dict(store.state["progress"] or dict(indexed=0,processed=0,comparisons=0,pending=None))
            def save():
                store.snapshot(db,progress)
            for ordinal,payload in enumerate(lines_from_paths(exact_paths),1):
                if ordinal <= progress["indexed"]:
                    continue
                completion=loads_completion(payload)
                db.execute("INSERT INTO exact VALUES (?,?,?)",(ordinal,repr(exact_completion_bucket_key(completion)),payload))
                progress["indexed"]=ordinal
                budget.boundary(save,progress=dict(stage="coverage-index",**progress))
            cache=DecodedCache()
            for ordinal,payload in enumerate(lines_from_paths(historical_paths),1):
                if ordinal <= progress["processed"]:
                    continue
                completion=loads_completion(payload)
                pending=progress["pending"]
                if pending is None:
                    try:
                        validate_completion(completion)
                    except InvalidMomentumGraph:
                        # Unknown or ambiguous physical routing is not a
                        # demonstrated vanishing historical class.
                        raise
                    except ValueError as error:
                        db.execute("INSERT INTO coverage VALUES (?,?,?)",(ordinal,"invalid",canonical_json(_invalid_historical_class(completion,error))))
                        progress["processed"]=ordinal
                        budget.boundary(save,progress=dict(stage="coverage",**progress))
                        continue
                    pending=dict(ordinal=ordinal,payload=payload,cursor=0)
                    progress["pending"]=pending
                elif pending["ordinal"] != ordinal or pending["payload"] != payload:
                    raise ValueError("Historical input changed during resume")
                match=None
                for sequence,known in db.execute("SELECT id,payload FROM exact WHERE bucket=? AND id>? ORDER BY id",
                    (repr(exact_completion_bucket_key(completion)),pending["cursor"])):
                    progress["comparisons"]+=1
                    pending["cursor"]=sequence
                    if are_equivalent_completions(completion,cache.get(known)):
                        match=sequence
                        break
                    budget.boundary(save,progress=dict(stage="coverage-comparison",processed=progress["processed"],comparisons=progress["comparisons"]))
                if match is not None:
                    status,evidence="reproduced",dict(historical_record=ordinal,generated_class=match)
                else:
                    audit=audit_amplitude_symmetrisation(completion)
                    if audit.status == "unsupported":
                        raise ValueError(f"UNSUPPORTED_HISTORICAL_AMPLITUDE: {audit.reason}")
                    if audit.is_zero:
                        status,evidence="invalid",_invalid_historical_class(completion,ValueError(audit.reason))
                    else:
                        status,evidence="missing",dict(code="MISSING_HISTORICAL_CLASS",record=payload)
                db.execute("INSERT INTO coverage VALUES (?,?,?)",(ordinal,status,canonical_json(evidence)))
                progress["processed"]=ordinal
                progress["pending"]=None
                budget.boundary(save,progress=dict(stage="coverage",processed=ordinal,comparisons=progress["comparisons"]))
            progress["counts"]=dict(db.execute("SELECT status,COUNT(*) FROM coverage GROUP BY status"))
            output=store.publish_lines(canonical_json(dict(record=i,status=status,evidence=__import__('json').loads(evidence)))
                for i,status,evidence in db.execute("SELECT id,status,evidence FROM coverage ORDER BY id"))
            store.save(progress,artifacts={"output":output})
            store.snapshot(db,progress,complete=True)
            return store.state
