"""Record-level old/current correspondence with explicit physical migration."""
import argparse
from collections import Counter,defaultdict
import json
from pathlib import Path

from neutrinomass.database.serialization import iter_completion_jsonl,dumps_completion,loads_completion
from neutrinomass.database.resumable_reduction import lines_from_paths
from neutrinomass.database.checkpoints import digest_file,atomic_json
from neutrinomass.database.rebuild import historical_record_normalizer,validate_completion,model_strings
from neutrinomass.completions.physical_routes import physical_numerators
from neutrinomass.completions.completions import _completion_with_momentum_contributions,exact_completion_bucket_key,are_equivalent_completions
from neutrinomass.completions.amplitudes import audit_amplitude_symmetrisation
from neutrinomass.database.global_equivalence import global_lagrangian_bucket,equivalent_global_lagrangians


def compare(operator,old_path,new_path):
    new=list(iter_completion_jsonl(new_path))
    buckets=defaultdict(list)
    for i,record in enumerate(new):
        buckets[exact_completion_bucket_key(record)].append((i,record))
    normalise=historical_record_normalizer(operator)
    rows=[];covered=set();old_interactions=defaultdict(list);old_models=defaultdict(list)
    for ordinal,payload in enumerate(lines_from_paths([old_path])):
        old=loads_completion(payload)
        old_interactions[global_lagrangian_bucket(old)].append((ordinal,old))
        old_models[tuple(model_strings(old))].append(ordinal)
        row=dict(old_record=ordinal,old_payload=payload,new_records=[])
        try:
            physical=physical_numerators(old)
            row["old_routes"]=repr(old.momentum_contributions)
            row["physical_routes"]=repr(physical)
            row["routing_changed"]=tuple(old.momentum_contributions)!=tuple(physical)
            if len(physical)!=1:
                row["code"]="INVALID_OLD_DERIVATIVE_BUDGET"
                row["physical_numerator_count"]=len(physical)
            else:
                migrated=_completion_with_momentum_contributions(old,physical)
                migrated.lorentz_projection=None
                migrated,=normalise(migrated)
                validate_completion(migrated,check_vanishing=False)
                zero=[i for i,t in enumerate(migrated.terms) if t.safe_simplify()==0]
                if zero:
                    row.update(code="DECODED_VANISHING_VERTEX",terms=zero)
                else:
                    audit=audit_amplitude_symmetrisation(migrated)
                    if audit.status=="unsupported":
                        raise ValueError("UNSUPPORTED_OLD_AMPLITUDE: "+audit.reason)
                    if audit.is_zero:
                        row["code"]="ZERO_PHYSICAL_AMPLITUDE"
                    else:
                        matches=[i for i,record in buckets[exact_completion_bucket_key(migrated)]
                                 if are_equivalent_completions(migrated,record)]
                        row.update(code="MATCHED_AFTER_PHYSICAL_NORMALIZATION" if matches else "UNEXPLAINED_OMISSION",
                                   new_records=matches,migrated_payload=dumps_completion(migrated))
                        covered.update(matches)
        except ValueError as error:
            row.update(code="BLOCKING_OLD_RECORD_DIAGNOSTIC",diagnostic=str(error))
        rows.append(row)
    additions=[dict(new_record=i,record=dumps_completion(record),model=model_strings(record),
                    old_interaction_matches=[ordinal for ordinal,old in old_interactions[global_lagrangian_bucket(record)]
                                             if equivalent_global_lagrangians(old,record)],
                    old_representation_matches=old_models[tuple(model_strings(record))])
               for i,record in enumerate(new) if i not in covered]
    return dict(schema_version=1,operator=operator,old_count=len(rows),new_count=len(new),
        old_sha256=digest_file(old_path),new_sha256=digest_file(new_path),
        counts=dict(Counter(row["code"] for row in rows)),correspondence=rows,added=additions,
        matched_current_classes=len(covered),
        duplicate_old_representatives=sum(bool(row["new_records"]) for row in rows)-len(covered),
        unexplained=sum(row["code"] in {"UNEXPLAINED_OMISSION","BLOCKING_OLD_RECORD_DIAGNOSTIC"} for row in rows))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operator");parser.add_argument("old",type=Path)
    parser.add_argument("new",type=Path);parser.add_argument("output",type=Path)
    args=parser.parse_args()
    result=compare(args.operator,args.old,args.new)
    atomic_json(args.output,result)
    print(json.dumps({k:v for k,v in result.items() if k not in {"correspondence","added"}},indent=2))
    raise SystemExit(bool(result["unexplained"]))
