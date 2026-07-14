#!/usr/bin/env bash

# Source this file from a Spartan login node. Configuration is via the
# NEUTRINOMASS_* variables documented in cluster/README.md.
priority4_submit_spartan() (
    set -euo pipefail

    local repo_root project_root account output_root legacy_dir task_manifest
    local source_commit python_cmd task_count array_job final_job
    repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
    project_root="$(cd "${repo_root}/.." && pwd)"
    account="${NEUTRINOMASS_ACCOUNT:-punim0011}"
    output_root="${NEUTRINOMASS_OUTPUT_ROOT:-${project_root}/priority4-rebuild-v5}"
    legacy_dir="${NEUTRINOMASS_LEGACY_DIR:-${project_root}/raw_completions}"
    task_manifest="${output_root}/cluster_tasks.json"
    source_commit="9d96add59fc636223eda1e39d2f16a2a5de5e908"

    if [[ "$(git -C "${repo_root}" branch --show-current)" != "cluster-rebuild" ]]; then
        echo "Check out branch cluster-rebuild before submission." >&2
        return 2
    fi
    if [[ -n "$(git -C "${repo_root}" status --porcelain --untracked-files=no)" ]]; then
        echo "Refusing to submit from a dirty tracked worktree." >&2
        return 2
    fi
    if ! command -v sbatch >/dev/null 2>&1; then
        echo "sbatch is unavailable; run this on a Spartan login node." >&2
        return 2
    fi

    source "${repo_root}/cluster/load_spartan.sh"
    python_cmd="${NEUTRINOMASS_PYTHON:-${repo_root}/.venv/bin/python}"
    if [[ ! -x "${python_cmd}" ]]; then
        if [[ -n "${NEUTRINOMASS_PYTHON:-}" ]]; then
            echo "Configured Python is not executable: ${python_cmd}" >&2
            return 2
        fi
        python3 -m venv "${repo_root}/.venv"
        "${python_cmd}" -m pip install --upgrade pip
        "${python_cmd}" -m pip install -r "${repo_root}/requirements-test.lock"
        "${python_cmd}" -m pip install --no-deps -e "${repo_root}"
    fi

    # Python 3.12 removed stdlib distutils, but the validated SymPy 1.2 stack
    # imports it. This exact setuptools version supplies the compatibility shim
    # used by the validated local Python 3.12 environment. Run this outside the
    # creation guard so an already-created Spartan .venv is repaired in place.
    "${python_cmd}" -m pip install "setuptools==69.5.1"

    mkdir -p "${output_root}/logs" "${output_root}/.matplotlib-submit"
    export MPLCONFIGDIR="${output_root}/.matplotlib-submit"
    cd "${repo_root}"
    "${python_cmd}" -m pytest -q \
        cluster_rebuild_test.py \
        cluster_legacy_archive_test.py \
        compare_published_database_test.py \
        rebuild_completion_database_test.py
    "${python_cmd}" cluster/fetch_legacy_archive.py "${legacy_dir}"

    "${python_cmd}" cluster_rebuild.py plan \
        "${legacy_dir}" \
        --source-commit "${source_commit}" \
        --output "${task_manifest}"
    task_count="$("${python_cmd}" cluster_rebuild.py count "${task_manifest}")"
    if [[ "${task_count}" != "486" ]]; then
        echo "Expected 486 tasks, found ${task_count}." >&2
        return 2
    fi

    array_job="$(sbatch \
        --parsable \
        --account "${account}" \
        --job-name nm-p4-census \
        --cpus-per-task 1 \
        --mem "${NEUTRINOMASS_MEMORY:-8G}" \
        --time "${NEUTRINOMASS_WALLTIME:-24:00:00}" \
        --array "0-$((task_count - 1))%${NEUTRINOMASS_CONCURRENCY:-16}" \
        --output "${output_root}/logs/census_%A_%a.out" \
        --error "${output_root}/logs/census_%A_%a.err" \
        --export "ALL,NEUTRINOMASS_REPO_ROOT=${repo_root},NEUTRINOMASS_TASK_MANIFEST=${task_manifest},NEUTRINOMASS_OUTPUT_ROOT=${output_root},NEUTRINOMASS_LEGACY_DIR=${legacy_dir}" \
        "${repo_root}/cluster/slurm_priority4_worker.sh")"
    array_job="${array_job%%;*}"

    final_job="$(sbatch \
        --parsable \
        --account "${account}" \
        --dependency "afterok:${array_job}" \
        --job-name nm-p4-finalize \
        --cpus-per-task 1 \
        --mem "${NEUTRINOMASS_FINAL_MEMORY:-16G}" \
        --time "${NEUTRINOMASS_FINAL_WALLTIME:-12:00:00}" \
        --output "${output_root}/logs/finalize_%j.out" \
        --error "${output_root}/logs/finalize_%j.err" \
        --export "ALL,NEUTRINOMASS_REPO_ROOT=${repo_root},NEUTRINOMASS_TASK_MANIFEST=${task_manifest},NEUTRINOMASS_OUTPUT_ROOT=${output_root},NEUTRINOMASS_LEGACY_DIR=${legacy_dir}" \
        "${repo_root}/cluster/slurm_priority4_finalize.sh")"
    final_job="${final_job%%;*}"

    printf 'Submitted Priority-4 census array %s and finalizer %s.\n' \
        "${array_job}" "${final_job}"
    printf 'Monitor with: squeue -j %s,%s\n' "${array_job}" "${final_job}"
)

if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
    priority4_submit_spartan "$@"
    return $?
fi
priority4_submit_spartan "$@"
