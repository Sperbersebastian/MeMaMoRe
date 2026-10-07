#!/usr/bin/env bash
# create_env_viruses.sh
# Umbrella creator for the virus module (bin/modules/viruses.sh).
# Own envs:    virsorter2_env, vibrant_env, checkv_env (per-tool creators)
# Shared envs: env_assembly_core, viralverify_env, genomad_env, plasmids_core
#              (provisioned by create_env_plasmids.sh, incl. the geNomad DB)
# Existing envs are kept; FORCE_CLEAN=1 rebuilds the virus-specific ones.
set -euo pipefail

export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
command -v micromamba >/dev/null 2>&1 || { echo "[err] micromamba not found"; exit 1; }

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FORCE_CLEAN="${FORCE_CLEAN:-0}"

env_exists(){ micromamba env list | awk 'NR>2{print $1}' | grep -qx "$1"; }

# --- shared envs from the plasmid module ---
for e in env_assembly_core viralverify_env genomad_env plasmids_core; do
  if ! env_exists "$e"; then
    echo "[deps] $e missing -> running create_env_plasmids.sh"
    bash "$ROOT_DIR/envs/create_env_plasmids.sh"
    break
  fi
done

# --- virus-specific envs ---
for tool in virsorter2 vibrant checkv; do
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

echo "[creator] ensured: env_assembly_core viralverify_env genomad_env plasmids_core virsorter2_env vibrant_env checkv_env"
