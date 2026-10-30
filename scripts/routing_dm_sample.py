"""Bounded package-generation check with distinct heavy edge species."""
from importlib import import_module
from itertools import islice
import json
from pathlib import Path
import sys

from neutrinomass.tensormethod import Field, L, H, D, eps
from neutrinomass.tensormethod.core import BOSE
from neutrinomass.completions import EffectiveOperator
from neutrinomass.completions.routing import validate_routes
from neutrinomass.database.serialization import dumps_completion

module = import_module("neutrinomass.completions.completions")
factory = module.exotic_field_and_term
module.exotic_field_and_term = lambda op, symbols, field_dict: factory(op, symbols, {})
phi = Field("phidm", "00000", comm=BOSE)
operator = EffectiveOperator("DM_scalar_n",
    L("u0 i0") * H("i1") * L.conj("d0 i2") * H.conj("i3")
    * eps("-i0 -i1") * eps("-i2 -i3") * phi("") * D(phi.conj, "11")("u1 d1"))
output = Path(sys.argv[1])
count = 0
with output.open("x") as handle:
    for completion in islice(module.momentum_routed_completion_stream(operator, canonical_partitions=True), 64):
        validate_routes(completion, derivative_order=1)
        handle.write(dumps_completion(completion) + "\n")
        count += 1
print(json.dumps({"sample_records": count, "complete_catalogue": False, "source": module.__file__}))
