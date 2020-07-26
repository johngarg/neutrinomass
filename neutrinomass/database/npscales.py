#!/usr/bin/env python3

from neutrinomass.database.utils import estimate_np_scale
from neutrinomass.completions import EFF_OPERATORS, DERIV_EFF_OPERATORS

for label, op in EFF_OPERATORS.items():
    print(f"{label}: {max(estimate_np_scale(op))}")
