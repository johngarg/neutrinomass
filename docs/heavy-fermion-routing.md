# Heavy-fermion routing development checkpoint

Branch `fix/heavy-fermion-routing`, based on
`d77831f645d926b3502bd71c76ad04050d7c20dd`. This is **not a release or a
certified replacement derivative census**.

`derivative_route_alternatives` now selects the undifferentiated heavy
endpoint transformed by `process_derivative_term`. For two heavy fermions,
the differentiated helper and physical numerator belong to different edges.
For one heavy and one light fermion, only the light helper is admissible.
The vectorlike transformation is partner conjugation, not conjugation alone.

Construction retains interaction indices and concrete heavy tensor ports on
graph nodes. Route resolution uses the input leaf occurrence, not the global
species dictionary. Equal occurrences that cannot be resolved from their
states and gauge indices raise an explicit unsupported-provenance exception.
The endpoint validator independently checks chirality, complete species
orbits, conjugate gauge states and vectorlike Weyl partners. It checks
numerator degree separately from denominator expansion order.

Routed construction is checked before historical projection. Consumers of
supplied route metadata (equivalence, completion fingerprints, amplitude
momentum expansion, JSON and executable export) check endpoints again.
There is no certification cache on mutable completions.

JSON version 3 preserves graph occurrence ports and validates supplied routes
on reading and writing. Versions 1 and 2 remain decodable for auditing;
decoding is not validation or migration. `migrate_legacy_one_derivative`
repairs only a complete, unambiguously reconstructed graph with one supplied
kinetic contribution, clears its old Lorentz projection, and leaves the input
unchanged. Keep original bytes, hashes and producer information separately.
It cannot recover candidates removed by earlier projection or deduplication.

## Evidence

The two original saved DM records are in `tests/fixtures/scalar-n-routing.jsonl`;
their source hash, line numbers and producer manifest are retained alongside
them. Record 115 has physical edge `(17496, 157464)` and zero scalar-current
projection at constant Higgs vev with massless external neutrinos. Record 116
remains nonzero. These are not full off-shell electroweak EFT identities.

The read-only audit script reproduces 441 discrepancies among 879 pruned
scalar records and 996 among 2,196 unpruned scalar records. The momentum-cut
test changes 120 nonzero to zero and 120 zero to nonzero in the pruned streams;
the corresponding unpruned counts are 288 and 288. This test does not include
the weak projector. All 492 fermion-DM controls across both scans agree.
The bounded generator produces 64 fresh endpoint-consistent scalar records.

The bounded D6a fixture changes from four routed records to two. Its two
wrong-cut records become equivalent to the two consistent records after
migration. The full existing D3 test changes 26 to 25 records and seven to
six exact classes: its second routed record is a wrong-cut duplicate of the
first. The local contributions themselves are unchanged in that comparison.
`scripts/routing_benchmarks.py` runs against either checkout to preserve the
before/after witnesses. These are implementation benchmark counts, not
publication-wide physical model counts.

## Unresolved scope

The legacy *local derivative* path can produce a physical momentum edge
without any route metadata. In the saved 745-class D20 structural artifact,
156 records have that omission; among the other 589, 360 routes disagree and
229 agree. Current consumer guards validate supplied routes, and do not yet
certify that local derivative path. Merely attaching a reconstructed edge
after `project_local` would not establish the correctness of its earlier
projection or candidate selection. That integration is required before this
branch can be considered a package-wide repair.

The stored D5b files inspected here contain 432 raw/43 reduced records, all
without their physical momentum-route metadata. They are not the requested
724-record coverage baseline. Do not equate these files with that benchmark.

The local reuse manifest inventories 290 derivative artifacts conservatively.
No production reuse is approved. Raw graph repair is allowed only as an
explicit intermediate stage; projection, equivalence, amplitude audits and
dependent filtering/comparisons must be redone. Any missing raw candidates or
uncertain earlier pruning requires regeneration. Remote coverage has not
been inspected. The DM field-content catalogue has not been regenerated in
this checkpoint, so no new claim is made that its 45 scalar contents remain
unchanged. Multi-derivative generation and all published counts remain subject
to the broader task's completion criteria.

Run with the pinned sibling interpreter and explicit source path:

```sh
PYTHONPATH=. ../neutrinomass/.venv311/bin/python -m pytest neutrinomass/completions/routing_test.py
PYTHONPATH=. ../neutrinomass/.venv311/bin/python scripts/routing_reuse_manifest.py
```

The vault report and machine-readable results are in
`.Data/exploding-ops/routing-20261030/` relative to the vault root. Test logs
record their own pass/failure outcomes; a development checkpoint must not be
described as a full passing scientific validation.
