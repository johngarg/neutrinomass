#!/usr/bin/env bash
set -euo pipefail

task_id="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
repo_root="${NEUTRINOMASS_REPO_ROOT:?NEUTRINOMASS_REPO_ROOT is required}"
task_manifest="${NEUTRINOMASS_TASK_MANIFEST:?NEUTRINOMASS_TASK_MANIFEST is required}"
output_root="${NEUTRINOMASS_OUTPUT_ROOT:?NEUTRINOMASS_OUTPUT_ROOT is required}"
legacy_dir="${NEUTRINOMASS_LEGACY_DIR:?NEUTRINOMASS_LEGACY_DIR is required}"

source "${repo_root}/cluster/load_spartan.sh"
python_cmd="${NEUTRINOMASS_PYTHON:-${repo_root}/.venv/bin/python}"
if [[ ! -x "${python_cmd}" ]]; then
    echo "Missing cluster environment: ${python_cmd}" >&2
    exit 2
fi

scratch_root="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
cd "${repo_root}"
"${python_cmd}" cluster_rebuild.py worker \
    "${task_manifest}" \
    "${task_id}" \
    "${output_root}" \
    "${legacy_dir}" \
    --scratch-root "${scratch_root}"
