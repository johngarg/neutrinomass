#!/usr/bin/env bash
set -euo pipefail

repo_root="${NEUTRINOMASS_REPO_ROOT:?NEUTRINOMASS_REPO_ROOT is required}"
task_manifest="${NEUTRINOMASS_TASK_MANIFEST:?NEUTRINOMASS_TASK_MANIFEST is required}"
legacy_dir="${NEUTRINOMASS_LEGACY_DIR:?NEUTRINOMASS_LEGACY_DIR is required}"
source_commit="${NEUTRINOMASS_SOURCE_COMMIT:?NEUTRINOMASS_SOURCE_COMMIT is required}"

source "${repo_root}/cluster/load_spartan.sh"
python_cmd="${NEUTRINOMASS_PYTHON:-${repo_root}/.venv/bin/python}"
if [[ ! -x "${python_cmd}" ]]; then
    echo "Missing cluster environment: ${python_cmd}" >&2
    exit 2
fi

scratch_root="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
export MPLCONFIGDIR="${scratch_root}/matplotlib-priority4-preflight"
mkdir -p "${MPLCONFIGDIR}"
cd "${repo_root}"

"${python_cmd}" -m pytest -q \
    cluster_rebuild_test.py \
    cluster_legacy_archive_test.py \
    compare_published_database_test.py \
    rebuild_completion_database_test.py
"${python_cmd}" cluster_rebuild.py plan \
    "${legacy_dir}" \
    --source-commit "${source_commit}" \
    --output "${task_manifest}"

task_count="$("${python_cmd}" cluster_rebuild.py count "${task_manifest}")"
if [[ "${task_count}" != "486" ]]; then
    echo "Expected 486 tasks, found ${task_count}." >&2
    exit 2
fi
source "${repo_root}/cluster/resource_tiers.sh"
normal_array="$("${python_cmd}" cluster_rebuild.py array "${task_manifest}" normal)"
high_array="$("${python_cmd}" cluster_rebuild.py array "${task_manifest}" high)"
if [[ "${normal_array}" != "${PRIORITY4_NORMAL_TASK_ARRAY}" ]]; then
    echo "Normal-resource task array differs from the pinned split." >&2
    exit 2
fi
if [[ "${high_array}" != "${PRIORITY4_HIGH_TASK_ARRAY}" ]]; then
    echo "High-resource task array differs from the pinned split." >&2
    exit 2
fi
printf 'Validated Priority-4 census manifest with %s tasks.\n' "${task_count}"
