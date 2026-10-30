"""Bounded historical-projection witnesses, usable with either checkout."""
import argparse
from importlib import import_module
import json
from pathlib import Path
import sys

from neutrinomass.completions.operators import DERIV_EFF_OPERATORS
from neutrinomass.database.serialization import completion_to_record

parser = argparse.ArgumentParser()
parser.add_argument("output", type=Path)
args = parser.parse_args()
module = import_module("neutrinomass.completions.completions")
original = module.partitions
module.partitions = lambda op, verbose=False: original(op, verbose=verbose)[144:145]
bounded = module.momentum_routed_completions(DERIV_EFF_OPERATORS["D6a"])
module.partitions = original
result = {"interpreter": sys.executable, "source": module.__file__,
          "D6a": [completion_to_record(c) for c in bounded]}
args.output.write_text(json.dumps(result, indent=2) + "\n")
print("D6a", len(bounded), [c.lorentz_projection.coordinates for c in bounded], flush=True)
all_d3 = module.deriv_operator_completions(DERIV_EFF_OPERATORS["D3"])
result["D3"] = [completion_to_record(c) for c in all_d3]
args.output.write_text(json.dumps(result, indent=2) + "\n")
print("D3", len(all_d3), flush=True)
