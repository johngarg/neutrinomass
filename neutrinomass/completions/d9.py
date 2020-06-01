#!/usr/bin/env python3

from neutrinomass.completions.operators import EFF_OPERATORS
from neutrinomass.completions.completions import (
    operator_completions,
    collect_completions,
    filter_completions,
)

completions = {5: [], 7: [], 9: []}

for label, op in EFF_OPERATORS.items():
    print(f"Working on {label}...")
    dim = op.mass_dimension
    if dim < 11:
        model = operator_completions(op)
        completions[dim] += model
    else:
        break

prev = {}
for dim, models in completions.items():
    collected = collect_completions(models)
    filtered = filter_completions(collected, prev)
    filtered_models[dim] = filtered
    prev = filtered
