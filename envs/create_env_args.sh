#!/usr/bin/env bash
# create_env_args.sh
# Umbrella creator for the ARG module (bin/modules/args.sh).
# Envs: deeparg_env (also used by summary_args), rgi_env, amrplusplus_env
# Existing envs are kept; FORCE_CLEAN=1 rebuilds them.
set -euo pipefail

export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
command -v micromamba >/dev/null 2>&1 || { echo "[err] micromamba not found"; exit 1; }

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FORCE_CLEAN="${FORCE_CLEAN:-0}"

env_exists(){ micromamba env list | awk 'NR>2{print $1}' | grep -qx "$1"; }

for tool in deeparg rgi amrplusplus; do
  env="${tool}_env"
  if [[ "$FORCE_CLEAN" == "1" ]] && env_exists "$env"; then
    micromamba env remove -y -n "$env"
  fi
  if env_exists "$env"; then
    echo "[ok] $env exists"
  else
    bash "$ROOT_DIR/envs/create_env_${tool}.sh"
  fi
done

echo "[creator] ensured: deeparg_env rgi_env amrplusplus_env"
