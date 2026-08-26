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
a small validation job, separate normal and high-resource census arrays, and a
finalizer behind `afterok` dependencies. The validation job runs the tests,
creates the 486-task manifest, and verifies the pinned resource split on a
compute node rather than using the Spartan login node for Python workloads.

An already prepared Python environment can be selected with
`NEUTRINOMASS_PYTHON=/absolute/path/to/python`. Set
`NEUTRINOMASS_SKIP_MODULES=1` only when that environment does not require a
Spartan Python module.

That one command creates `.venv` when needed, downloads and verifies the
published 199 MB Zenodo `raw_completions.zip` when the legacy directory is
absent, safely extracts its 243 inputs (about 4.4 GB), submits a validation job
that writes a 486-task manifest (243 operators times two hash seeds), submits
two throttled Slurm arrays dependent on that validation, and submits an
`afterok` finalizer. Each array job uses node-local scratch and gzip-packages
the large raw, structural-exact, and physical-exact JSONL files before
publishing only its own operator directory. The
finalizer validates every checksum and historical comparison, compares both
seed reports, writes `migration_manifest.json`, regenerates the
democratic/one-loop-Weinberg filtering artifacts, and writes
`published_comparison.{json,csv,md}` with global and operator-level changes
from the counts reported in 2009.13537.

The normal tier defaults to one CPU, 8 GiB, 24 hours, and at most 16
simultaneous jobs. The 15 operators that exhausted 8 GiB in the completed v7
diagnostic run are submitted for both hash seeds in a separate tier with one
CPU, 16 GiB, 48 hours, and at most four simultaneous jobs. The validation job
defaults to one CPU, 8 GiB, and two hours. Override them before sourcing when
required:

```bash
export NEUTRINOMASS_ACCOUNT=punim0011
export NEUTRINOMASS_CONCURRENCY=12
export NEUTRINOMASS_MEMORY=8G
export NEUTRINOMASS_WALLTIME=24:00:00
export NEUTRINOMASS_HIGH_CONCURRENCY=4
export NEUTRINOMASS_HIGH_MEMORY=16G
export NEUTRINOMASS_HIGH_WALLTIME=48:00:00
export NEUTRINOMASS_LEGACY_DIR=/data/gpfs/projects/punim0011/garj/exploding-ops/raw_completions
export NEUTRINOMASS_OUTPUT_ROOT=/data/gpfs/projects/punim0011/garj/exploding-ops/priority4-rebuild-v10
source cluster/submit_spartan.sh
```

The high-resource operators are `71p`, `77p`, `78p`, `79a`, `79b`, `7p`,
`80a`--`80d`, `81a`--`81d`, and `8pp`. Earlier runs of `71p` and `79a`
reached the 16 GiB allocation and the 48-hour limit while materialising and
re-auditing tens of thousands of historical records. The historical census is
now streamed through a disk-backed exact-class index, and an already matched
class inherits the amplitude result of its audited exact survivor. V9 showed a
second list-backed bottleneck in raw regular-operator partition generation:
seven tasks again reached essentially 16 GiB and 48 hours. That partition
stream is now consumed one furnished graph at a time. Census phase and record
progress is also forwarded to the persistent Slurm output. The 16 GiB,
48-hour tier remains as headroom for these large inputs.

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

The workflow is restartable: sourcing the submitter again creates new arrays,
but valid completed task reports are checksum-verified and skipped. The source
commit for the scientific census is pinned to
`ab84e6d3d279a562c42e008b862c3416e1040613`. This merge pins scientific
commit `d77831f645d926b3502bd71c76ad04050d7c20dd`. In addition to the earlier
propagator, tensor-statistics, and amplitude checks, it streams and
disk-deduplicates legacy records, avoids redundant historical amplitude
audits, streams raw regular-operator partitions, and restricts pre-furnishing
derivative-partition canonicalisation to operators checked against a complete
raw baseline. The latter restores the missing D5b class while retaining the
verified D20 optimisation. Later commits on the cluster branch may add only
orchestration around those scientific sources.
