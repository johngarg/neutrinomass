#!/usr/bin/env bash
set -euo pipefail

repo_root="${NEUTRINOMASS_REPO_ROOT:?NEUTRINOMASS_REPO_ROOT is required}"
task_manifest="${NEUTRINOMASS_TASK_MANIFEST:?NEUTRINOMASS_TASK_MANIFEST is required}"
output_root="${NEUTRINOMASS_OUTPUT_ROOT:?NEUTRINOMASS_OUTPUT_ROOT is required}"
legacy_dir="${NEUTRINOMASS_LEGACY_DIR:?NEUTRINOMASS_LEGACY_DIR is required}"

source "${repo_root}/cluster/load_spartan.sh"
python_cmd="${repo_root}/.venv/bin/python"
manifest="${output_root}/migration_manifest.json"
filtered_dir="${output_root}/filtered"

cd "${repo_root}"
"${python_cmd}" cluster_rebuild.py consolidate \
    "${task_manifest}" \
    "${output_root}" \
    "${legacy_dir}" \
    --output "${manifest}"
"${python_cmd}" rebuild_completion_database.py filter \
    "${output_root}" \
    "${manifest}" \
    --output-dir "${filtered_dir}"
