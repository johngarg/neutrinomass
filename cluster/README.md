# Priority-4 Spartan rebuild

On a Spartan login node, check out `cluster-rebuild` and run:

```bash
source cluster/submit_spartan.sh
```

To prepare the environment and verified legacy inputs, but stop before any
Slurm submission, run:

```bash
export NEUTRINOMASS_PREPARE_ONLY=1
source cluster/submit_spartan.sh
unset NEUTRINOMASS_PREPARE_ONLY
```

After reviewing the prepared inputs, source the submitter normally. It submits
a small validation job, then the census array and finalizer behind `afterok`
dependencies. The validation job runs the tests and creates the 486-task
manifest on a compute node rather than using the Spartan login node for Python
workloads.

An already prepared Python environment can be selected with
`NEUTRINOMASS_PYTHON=/absolute/path/to/python`. Set
`NEUTRINOMASS_SKIP_MODULES=1` only when that environment does not require a
Spartan Python module.

That one command creates `.venv` when needed, downloads and verifies the
published 199 MB Zenodo `raw_completions.zip` when the legacy directory is
absent, safely extracts its 243 inputs (about 4.4 GB), submits a validation job
that writes a 486-task manifest (243 operators times two hash seeds), submits a
throttled Slurm array dependent on that validation, and submits an `afterok`
finalizer. Each array job uses node-local scratch and gzip-packages the large
raw, structural-exact, and physical-exact JSONL files before publishing only
its own operator directory. The
finalizer validates every checksum and historical comparison, compares both
seed reports, writes `migration_manifest.json`, regenerates the
democratic/one-loop-Weinberg filtering artifacts, and writes
`published_comparison.{json,csv,md}` with global and operator-level changes
from the counts reported in 2009.13537.

Defaults are suitable for the `punim0011` Spartan project: one CPU, 8 GiB,
24 hours, and at most 16 simultaneous census jobs. The validation job defaults
to one CPU, 8 GiB, and two hours. Override them before sourcing when required:

```bash
export NEUTRINOMASS_ACCOUNT=punim0011
export NEUTRINOMASS_CONCURRENCY=12
export NEUTRINOMASS_MEMORY=8G
export NEUTRINOMASS_WALLTIME=24:00:00
export NEUTRINOMASS_LEGACY_DIR=/data/gpfs/projects/punim0011/garj/exploding-operators/raw_completions
export NEUTRINOMASS_OUTPUT_ROOT=/data/gpfs/projects/punim0011/garj/exploding-operators/priority4-rebuild-v5
source cluster/submit_spartan.sh
```

The download is pinned to Zenodo record `4054618`, byte size `198689883`, and
published MD5 `f7a199f7718607e3740137e85c1d488b`. A cached archive is reused only
after both checks pass. An existing complete 243-file directory is reused; a
partial directory is never overwritten automatically.

The default Spartan module stack is `GCC/13.3.0 OpenBLAS/0.3.27
Python/3.12.3`; Python must be loaded after its compiler dependencies. Override
the complete ordered list with `NEUTRINOMASS_MODULES` if Spartan changes the
available versions.

The bootstrap explicitly installs `setuptools==69.5.1`. Python 3.12 no longer
ships `distutils`, but pinned SymPy 1.2 imports it; this setuptools release
provides the same compatibility shim as the validated local environment. The
install is deliberately repeated on submission so a partially created `.venv`
is repaired without manual deletion.

The workflow is restartable: sourcing the submitter again creates a new array,
but valid completed task reports are checksum-verified and skipped. The source
commit for the scientific census is pinned to
`9d96add59fc636223eda1e39d2f16a2a5de5e908`. This merge pins scientific
commit `4cd174efc1138c607d50f2848c5781dff6ebb878` and includes the repaired
propagator expansion, the unreduced no-EOM second-derivative placement bases,
streamed derivative generation, decoded-record vertex validation, and the
amplitude-level identical-field symmetrisation audit. Later commits on the
cluster branch may add only orchestration around those scientific sources.
