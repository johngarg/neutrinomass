"""Measure bounded-cache costs on identical immutable completion inputs.

Run each mode in a fresh interpreter under both seeds. Disabling memoization
changes no equivalence test, input order or representative survival ranking.
Compare the output artifact SHA-256 and ordered fingerprints across modes.
"""
import argparse
import json
from pathlib import Path
import time

from neutrinomass.completions import equivalence
from neutrinomass.database import resumable_reduction
from neutrinomass.database.checkpoints import WorkBudget,atomic_json,digest_file,memory_sample
from neutrinomass.database.global_equivalence import DecodedCache
from neutrinomass.database.provenance import runtime_provenance,scientific_contract
from neutrinomass.database.rebuild import amplitude_survival_rank


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input",type=Path);parser.add_argument("output",type=Path)
    parser.add_argument("--mode",choices=("bounded","disabled"),required=True)
    parser.add_argument("--derivative-ranking",action="store_true")
    args=parser.parse_args()
    provenance=runtime_provenance()
    if args.mode=="disabled":
        resumable_reduction.DecodedCache=lambda:DecodedCache(max_records=0)
        equivalence._GRAPH_CACHE_MAXSIZE=0
    equivalence.clear_interaction_graph_cache()
    if args.output.exists():raise ValueError("Benchmark output must be fresh")
    options=dict(representative_rank=amplitude_survival_rank,maximum_rank=1) if args.derivative_ranking else {}
    started=time.monotonic();cpu=time.process_time()
    result=resumable_reduction.reduce_resumable([args.input],args.output,
        dict(scientific_contract(provenance),benchmark_cache_mode=args.mode),budget=WorkBudget(),**options)
    summary=dict(producer=provenance,validator_sha256=digest_file(__file__),cache_mode=args.mode,
        input=dict(path=str(args.input),sha256=digest_file(args.input)),
        elapsed_seconds=time.monotonic()-started,cpu_seconds=time.process_time()-cpu,
        memory=memory_sample(),progress=result["progress"],output=result["artifacts"]["output"])
    atomic_json(args.output/"benchmark.json",summary)
    print(json.dumps(summary,sort_keys=True))


if __name__=="__main__":main()
