# Priority-4 Spartan rebuild

On a Spartan login node, check out `cluster-rebuild`, place the preserved
legacy `raw_completions/` directory beside the repository, and run:

```bash
source cluster/submit_spartan.sh
```

That one command creates `.venv` when needed, verifies the 243-file legacy
inventory, writes a 486-task manifest (243 operators times two hash seeds),
submits a throttled Slurm array, and submits an `afterok` finalizer. Each array
job uses node-local scratch, gzip-packages the large raw and exact JSONL files,
and publishes only its own operator directory. The finalizer validates every
checksum and historical comparison, compares both seed reports, writes
`migration_manifest.json`, and regenerates the democratic/one-loop-Weinberg
filtering artifacts.

Defaults are suitable for the `punim0011` Spartan project: one CPU, 8 GiB,
24 hours, and at most 24 simultaneous census jobs. Override them before
sourcing when required:

```bash
export NEUTRINOMASS_ACCOUNT=punim0011
export NEUTRINOMASS_CONCURRENCY=12
export NEUTRINOMASS_MEMORY=8G
export NEUTRINOMASS_WALLTIME=24:00:00
export NEUTRINOMASS_LEGACY_DIR=/data/gpfs/projects/punim0011/garj/exploding-operators/raw_completions
export NEUTRINOMASS_OUTPUT_ROOT=/data/gpfs/projects/punim0011/garj/exploding-operators/priority4-rebuild-v4
source cluster/submit_spartan.sh
```

The workflow is restartable: sourcing the submitter again creates a new array,
but valid completed task reports are checksum-verified and skipped. The source
commit for the scientific census is pinned to
`ce3caefd8fdba1535f9b7e373423d9c933e057b9`; the cluster branch may add only
orchestration around those unchanged scientific sources.
