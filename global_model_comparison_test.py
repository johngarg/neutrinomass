import json
import pytest
from scripts.compare_global_models import compare,sources_from_report,validate_totals
from neutrinomass.database.checkpoints import WorkBudget,BudgetExpired,digest_file
from neutrinomass.database.resumable_reduction import lines_from_paths


def test_global_model_union_recovers_without_counting_operator_transfers_as_additions(tmp_path):
    first,second="S,00,0,7,0","S,00,0,8,0"
    sources=[]
    def add(name,mode,side,rows,operator=None,legacy=False):
        path=tmp_path/name
        path.write_text("\n".join(rows)+"\n")
        sources.append(dict(mode=mode,side=side,operator=operator,kind="legacy" if legacy else "json",
            artifact=dict(path=str(path),sha256=digest_file(path))))
    add("old-a", "unfiltered","historical",[f"LazyCompletion(head={dict(operator_name='a',quantum_numbers=[first])!r}, tail='unused')"],"a",True)
    # This model moves from operator a to b and also occurs twice in b.
    add("new-b","unfiltered","current",[json.dumps(dict(fields=[field])) for field in (first,first,second)],"b")
    add("old-filter","published","historical",[json.dumps(dict(operator="a",fields=[first]))])
    add("new-filter","published","current",[json.dumps(dict(operator="b",fields=[first]))])
    full=compare(sources,tmp_path/"full",{})
    with pytest.raises(BudgetExpired):compare(sources,tmp_path/"resume",{},budget=WorkBudget(fault_after=2))
    resumed=compare(sources,tmp_path/"resume",{})
    assert resumed["progress"]==full["progress"]
    assert resumed["progress"]["modes"]["unfiltered"]==dict(historical=1,current=2,unchanged=1,added=1,removed=0)
    assert resumed["progress"]["modes"]["published"]==dict(historical=1,current=1,unchanged=1,added=0,removed=0)
    output=tmp_path/"resume"/resumed["artifacts"]["output"]["path"]
    rows=list(map(json.loads,lines_from_paths([output])))
    assert len(rows)==1 and rows[0]["fields"]==[second] and rows[0]["code"]=="MODEL_ADDED"


def test_global_comparison_does_not_publish_partial_inventory_totals(tmp_path):
    report=tmp_path/"report.json";report.write_text(json.dumps(dict(complete_census=False)))
    manifest=tmp_path/"manifest.json";manifest.write_text("{}")
    with pytest.raises(ValueError,match="COMPLETE_CENSUS"):
        sources_from_report(report,manifest)


def test_global_comparison_requires_agreement_with_independent_totals():
    baseline=dict(conjugacy_unfiltered_models=141989,filtered_models=11216)
    report=dict(global_representation_sets=142000,filtered_published=dict(representation_sets=11000))
    counts=dict(unfiltered=dict(historical=141989,current=142000),
                published=dict(historical=11216,current=11000))
    validate_totals(counts,report,baseline)
    counts["published"]["historical"]+=15
    with pytest.raises(ValueError,match="TOTAL_MISMATCH"):
        validate_totals(counts,report,baseline)
