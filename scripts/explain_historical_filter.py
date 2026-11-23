"""Save direct loop witnesses for the 11,231 versus 11,216 archive discrepancy."""
import argparse
from itertools import combinations
import json
from pathlib import Path

from neutrinomass.database.checkpoints import atomic_json,digest_file
from neutrinomass.database.filtering import published_model_key,published_one_loop_templates
from neutrinomass.database.provenance import runtime_provenance
from neutrinomass.database.rebuild import generated_one_loop_weinberg_models


def explain(baseline,output):
    baseline=Path(baseline);data=json.loads(baseline.read_text())
    if not data["passed"] or data["per_operator_differences"]:
        raise ValueError("HISTORICAL_BASELINE_NOT_VERIFIED")
    published=dict(published_one_loop_templates())
    v10=generated_one_loop_weinberg_models()
    literal=set(v10);conjugate={published_model_key(fields):fields for fields in v10}
    rows=[];distinct=set()
    for operator,delta in sorted(data["v10_filter_differences"].items()):
        if delta["published_only"]:raise ValueError("UNEXPLAINED_HISTORICAL_REMOVAL")
        for fields in delta["v10_only"]:
            key=published_model_key(fields)
            subsets=[s for n in range(1,len(key)+1) for s in combinations(key,n)]
            old_match=next((s for s in subsets if s in published),None)
            new_match=next((s for s in subsets if s in conjugate),None)
            raw_match=any(s in literal for n in range(1,len(fields)+1) for s in combinations(sorted(fields),n))
            if old_match is None or new_match is None or raw_match:
                raise ValueError("FILTER_DISCREPANCY_NEEDS_ADDITIONAL_EXPLANATION")
            distinct.add(tuple(fields))
            rows.append(dict(operator=operator,fields=fields,code="MISSED_CONJUGATE_LOOP_IDENTITY",
                published_loop=published[old_match],v10_loop_before_conjugacy=conjugate[new_match],
                published_loop_key=old_match,canonical_v10_loop_key=new_match,
                v10_literal_match=False,v10_conjugacy_match=True))
    if len(distinct)!=data["v10_filtered_models"]-data["filtered_models"]:
        raise ValueError("FILTER_DISCREPANCY_GLOBAL_COUNT_DOES_NOT_CLOSE")
    result=dict(passed=True,schema_version=1,environment=runtime_provenance(),
        validator_sha256=digest_file(__file__),baseline=dict(path=str(baseline.resolve()),sha256=digest_file(baseline)),
        operator_incidences=len(rows),distinct_added_models=len(distinct),
        per_operator={op:len(delta["v10_only"]) for op,delta in data["v10_filter_differences"].items()},
        explanation="Every extra v10 model matches a generated v10 loop template after conjugate identification, "
                    "but none matches a literal v10 template subset. Archived baryon lookup and repeated template "
                    "factors are separate historical-contract requirements.",witnesses=rows)
    atomic_json(output,result)
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline",type=Path);parser.add_argument("output",type=Path)
    args=parser.parse_args();result=explain(args.baseline,args.output)
    print(json.dumps({k:v for k,v in result.items() if k not in {"environment","witnesses"}},sort_keys=True))
