#!/usr/bin/env bash

# This file is sourced by both the login-shell submitter and Slurm workers.
if [[ "${NEUTRINOMASS_SKIP_MODULES:-0}" == "1" ]]; then
    # shellcheck disable=SC2317
    return 0 2>/dev/null || exit 0
fi

if ! command -v module >/dev/null 2>&1; then
    echo "Spartan environment-modules command is unavailable." >&2
    # This file is normally sourced, but remains safe when executed directly.
    # shellcheck disable=SC2317
    return 1 2>/dev/null || exit 1
fi

read -r -a priority4_modules <<< "${NEUTRINOMASS_MODULES:-Python/3.12.3}"
for priority4_module in "${priority4_modules[@]}"; do
    if ! module is-loaded "${priority4_module}" >/dev/null 2>&1; then
        module load "${priority4_module}"
    fi
done
unset priority4_module priority4_modules

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
