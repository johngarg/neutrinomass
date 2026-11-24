# Reliable erratum census

This workflow covers the published 243 scalar/fermion operators under hash
seeds 0 and 1. It excludes the separate vector extension. The full 71p pilot
and complete corrected census are still in progress as of 23 November 2026;
passing a regression suite is not a completed census.

## Scientific contracts

The reports keep five different quantities separate:

| Quantity | Identification |
| --- | --- |
| Generated records | Ordered output of the recorded generation ranges |
| Physical completion classes | Exact completion equivalence after decoded vertex, physical route and symmetrised amplitude validation |
| Global interaction Lagrangians | Exact interaction-multiset equivalence under species, conjugate and Dirac naming conventions; contractions and multiplicities retained; operator, topology and EFT routing excluded |
| Representation sets | The paper's democratic identification convention |
| Species multisets | Actual admitted copies, particle nature and baryon assignments |

Fingerprints prioritize exact comparisons; they do not prove equality.
Required derivative records must have explicit physical propagator routing.
Kinetic numerators follow the UV endpoints, and remaining even derivative
order is distributed over propagator denominators. Unsupported required
physics blocks publication.

`published` filtering preserves the matched historical convention.
`species-aware` filtering searches additional trees and heavy T1/T3 Weinberg
loops using only the candidate's admitted copies, with species assigned
before vanishing-vertex rejection. Exclusions carry concrete witnesses.
Both assume generic allowed couplings, the candidate's baryon assignments
and no additional imposed symmetries. Neither proves the stronger low-scale
genuineness claim. Version-two reports carry explicit TeV units.

The reproduced historical baseline is 141,989 unfiltered and 11,216 filtered
democratic models, including matching per-operator filtered sets. A faithful
replay of the old contraction-blind interaction key gives 430,810 and 11,483
under its respective unfiltered and first-representative filtered counting
rules. These are historical baselines, not corrected census totals.
`scripts/explain_historical_filter.py` explains the separate 11,231 result:
15 global models, appearing in 40 entries of 27b/27d/29a/29b, had conjugate
one-loop counterparts missed by the literal v10 lookup.

## Stages and continuations

`census_operator.py` retains whole-operator invocation and uses the same
staged implementation for generation, decoded validation, local/structural
reduction, amplitudes, historical normalization/reduction, coverage and
packaging. The persistent plan records stable topology/assignment ranges.
Use `--unit N` for generation units and `--stage-unit N` for amplitude or
coverage shards; `--topology`, `--producer`, `--placement`, `--colour` and
aligned assignment bounds can select generation ranges explicitly.

Inside a one-CPU `sapphire` Slurm allocation under `punim0011`, for example:

```bash
module purge
module load GCCcore/11.3.0 Python/3.11.3
ERRATUM_RUN=/data/gpfs/projects/punim0011/garj/exploding-ops/erratum-20261123
cd "$ERRATUM_RUN/worktrees/451530f"
export PYTHONPATH="$PWD" PYTHONHASHSEED=0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$ERRATUM_RUN/venv/bin/python" census_operator.py 71p \
  "$ERRATUM_RUN/historical/op_71p.dat" \
  "$ERRATUM_RUN/pilot-451530f/seed-0/71p" \
  --stage structural --resume --time-budget 3300 \
  --scratch "$SLURM_TMPDIR" --require-slurm
```

The active queue owns its units: do not launch this manual command against
an owned stage. Exit 75 means a durable continuation, not success. Resume
with identical inputs and source/environment contracts. Incompatible or
corrupt checkpoints require diagnosis; temporary unpublished files are not
accepted as results. Completed stages are immutable. Persistent compressed
chunks and atomic manifests include consistent SQLite snapshots, including
the representative cursor inside an unfinished exact-reduction bucket.

## Gated campaign

The November run freezes scientific generation at `451530f`; the external
worker driver at `c281b64` repairs shared-plan startup without changing that
scientific package. Finalisation uses `6909818`; the chain/recovery scripts
use `285fe12`. The global model comparison is at `c3c22f8`. Run configurations,
installed dependency pins, lockfile hashes and imported package paths are
recorded independently of the moving development branch.

GitHub Actions uses Python 3.11.3 and the same dependency lock as Spartan.
The inherited Python 3.12 workflow rejected archived LaTeX backslash escapes
under pytest's warnings-as-errors policy. Supporting that newer interpreter
requires explicit legacy-reader migration; the erratum campaign does not
claim Python 3.12 compatibility. The historical inputs and warning checks
remain unchanged.

`scripts/run_erratum_campaign.py` follows the pilot, real process-loss
replicas, `scripts/validate_71p_gate.py`, execution/reuse manifests, the
remaining census and finalisation. The gate requires both complete
264-topology 71p seeds, the 758,208 raw-generation reference, exact scientific
agreement, complete historical accounting, verified recovery/substitution
and measured resource headroom. No old completion stage is reusable merely
because an old report exists.

The authoritative campaign ledger includes completed allocations and the
maximum allocation of submitted jobs. It enforces 3,000 core-hours overall,
200 initially for validation/pilot, four workers, one CPU and at most two
hours per worker. Normal work units target 55 minutes; workers begin their
graceful exit before 90 minutes. Only diagnosed infrastructure failures
allow up to two automatic retries. Scientific failures, corruption and
repeated continuations without progress block submission.

After the gate, `scripts/prepare_erratum_campaign.py` records the exact reuse
and recomputation inventory. Finalisation is independently runnable through
`scripts/finalise_erratum.py` and `erratum_finalise_worker.py`; it does not
depend on an obsolete `afterok` chain. It refuses complete totals without
verified coverage of all 243 operators under both seeds.

`scripts/compare_global_models.py REPORT INPUT_MANIFEST OUTPUT` is a
recoverable, budgeted final follow-up. It compares global unions under the
published convention, so a model moving between operators is not called a
global addition/removal. Its totals must agree with the independent baseline
and final report. Its outputs supplement the per-operator changes, exact
global Lagrangian artifacts and filtering witnesses.

The live run's manifests, batch scripts, source changes, validation evidence
and ledger live beneath `erratum-20261123/`. The vault report
`20261123T151113.md` records measured evidence and remaining limitations.
Public release and erratum submission are outside this workflow.
