"""Stable class ranges and verified ordered publication of late-stage shards."""
from collections import Counter
from pathlib import Path

from neutrinomass.database.checkpoints import CheckpointStore,WorkBudget
from neutrinomass.database.resumable_reduction import input_contract,lines_from_paths
from neutrinomass.database.stages import map_records
from neutrinomass.database.coverage import historical_coverage,build_coverage_index


def ranges(count,size):
    if count<0 or size<1:
        raise ValueError("Invalid stage sharding bounds")
    return [(start,min(start+size,count)) for start in range(0,count,size)] or [(0,0)]


def _path(directory,state,key="output"):
    return Path(directory)/state["artifacts"][key]["path"]


def split_records(paths,directory,contract,*,size=128,budget=None):
    """Read a long stream once into immutable ordered class chunks."""
    budget=budget or WorkBudget()
    identity=dict(contract,algorithm="record-shards-v1",size=size,inputs=input_contract(paths))
    with CheckpointStore(directory,identity) as store:
        if store.state["complete"]:
            return store.state
        progress=dict(store.state["progress"] or dict(processed=0,chunks=0))
        buffered=[]
        def publish():
            artifact=store.publish_lines(buffered)
            key=f"chunk-{progress['chunks']:08d}"
            progress["processed"]+=len(buffered)
            progress["chunks"]+=1
            store.save(progress,artifacts={key:artifact})
            buffered.clear()
            budget.boundary(lambda:store.save(progress),progress=dict(stage="record-shards",**progress))
        for ordinal,line in enumerate(lines_from_paths(paths),1):
            if ordinal<=progress["processed"]:
                continue
            buffered.append(line)
            if len(buffered)==size:
                publish()
        if buffered or not progress["chunks"]:
            publish()
        store.save(progress,complete=True)
        return store.state


def sharded_map(paths,count,directory,contract,transform,*,algorithm,budget=None,scratch=None,
                shard_size=128,selected_unit=None):
    directory=Path(directory)
    bounds=ranges(count,shard_size)
    if selected_unit is not None and selected_unit not in range(len(bounds)):
        raise ValueError("Invalid amplitude stage unit")
    split=split_records(paths,directory/"inputs",contract,size=shard_size,budget=budget)
    if split["progress"]["processed"]!=count:
        raise ValueError("AMPLITUDE_INPUT_COUNT_MISMATCH")
    states=[];outputs=[];evidence=[]
    for i,(start,stop) in enumerate(bounds):
        if selected_unit is not None and i!=selected_unit:
            continue
        base=directory/"shards"/f"{i:08d}"
        chunk=_path(directory/"inputs",split,f"chunk-{i:08d}")
        state=map_records([chunk],base,contract,transform,algorithm=algorithm,budget=budget,scratch=scratch,
                          record_start=start,record_stop=stop,ordinal_offset=start)
        if selected_unit is not None:
            return state
        states.append((start,state))
        outputs.append(_path(base,state));evidence.append(_path(base,state,"evidence"))
    identity=dict(contract,algorithm=algorithm+"-shard-merge-v1",ranges=bounds,
                  inputs=input_contract(outputs+evidence))
    with CheckpointStore(directory,identity) as store:
        if store.state["complete"]:
            return store.state
        progress=dict(processed=sum(s["progress"]["processed"]-start for start,s in states),
            accepted=sum(s["progress"]["accepted"] for _,s in states),
            rejected=sum(s["progress"]["rejected"] for _,s in states))
        if progress["processed"]!=count:
            raise ValueError("INCOMPLETE_AMPLITUDE_SHARD_COVERAGE")
        store.save(progress,artifacts={"output":store.publish_lines(lines_from_paths(outputs)),
            "evidence":store.publish_lines(lines_from_paths(evidence))},complete=True)
        return store.state


def sharded_coverage(exact_paths,historical_paths,count,directory,contract,*,budget=None,scratch=None,
                     shard_size=128,selected_unit=None):
    directory=Path(directory)
    bounds=ranges(count,shard_size)
    if selected_unit is not None and selected_unit not in range(len(bounds)):
        raise ValueError("Invalid coverage stage unit")
    index=build_coverage_index(exact_paths,directory/"index",contract,budget=budget,scratch=scratch)
    index_path=_path(directory/"index",index,"database")
    split=split_records(historical_paths,directory/"inputs",contract,size=shard_size,budget=budget)
    if split["progress"]["processed"]!=count:
        raise ValueError("HISTORICAL_INPUT_COUNT_MISMATCH")
    outputs=[];counts=Counter();processed=0;comparisons=0
    for i,(start,stop) in enumerate(bounds):
        if selected_unit is not None and i!=selected_unit:
            continue
        base=directory/"shards"/f"{i:08d}"
        chunk=_path(directory/"inputs",split,f"chunk-{i:08d}")
        state=historical_coverage(exact_paths,[chunk],base,contract,budget=budget,scratch=scratch,
                                 record_start=start,record_stop=stop,index_path=index_path,ordinal_offset=start)
        if selected_unit is not None:
            return state
        processed+=state["progress"]["processed"]-start
        comparisons+=state["progress"]["comparisons"]
        counts.update(state["progress"]["counts"])
        outputs.append(_path(base,state))
    identity=dict(contract,algorithm="coverage-shard-merge-v1",ranges=bounds,inputs=input_contract(outputs))
    with CheckpointStore(directory,identity) as store:
        if store.state["complete"]:
            return store.state
        if processed!=count:
            raise ValueError("INCOMPLETE_HISTORICAL_SHARD_COVERAGE")
        store.save(dict(processed=processed,comparisons=comparisons,counts=dict(counts)),
            artifacts={"output":store.publish_lines(lines_from_paths(outputs))},complete=True)
        return store.state
