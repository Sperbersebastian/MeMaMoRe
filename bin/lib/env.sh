#!/usr/bin/env bash
set -euo pipefail

# Micromamba init (idempotent)
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
eval "$(micromamba shell hook --shell=bash)"

_env_exists() {
  micromamba env list | awk '{print $1}' | grep -qx "$1"
}

_recreate_env() {
  local env_name="$1" create_script="$2"
  echo "[env] removing $env_name"
  micromamba env remove -n "$env_name" -y || true
  echo "[env] creating $env_name via $create_script"
  bash "$create_script"
}

# Force controls (exported by main.sh)
FORCE_ENV_MODE="${FORCE_ENV_MODE:-none}"   # none|all|one
FORCE_ENV_TARGET="${FORCE_ENV_TARGET:-}"   # env name or module key

# Map: module_key -> "<env_name> <create_script_rel_to_ROOT>"
_env_spec() {
  case "$1" in
    ingest)    echo "env_ingest    envs/create_env_ingest.sh" ;;
    qc)        echo "env_qc        envs/create_env_qc.sh" ;;
    assembly)  echo "env_assembly  envs/create_env_assembly.sh" ;;
    binning)   echo "env_binning   envs/create_env_binning.sh" ;;
    sra_tools) echo "env_sra_tools envs/create_env_sra_tools.sh" ;; # used by merge()
    *)         echo "UNKNOWN UNKNOWN"; return 1 ;;
  esac
}

ensure_env() {
  local env_name="$1" create_script="$2"
  local force_this="no"
  [[ "$FORCE_ENV_MODE" == "all" ]] && force_this="yes"
  if [[ "$FORCE_ENV_MODE" == "one" && "${FORCE_ENV_TARGET:-}" == "$env_name" ]]; then
    force_this="yes"
  fi

  if [[ "$force_this" == "yes" ]]; then
    _recreate_env "$env_name" "$create_script"
  elif ! _env_exists "$env_name"; then
    echo "[env] creating $env_name via $create_script"
    bash "$create_script"
  fi

  _env_exists "$env_name" || { echo "[env] failed to create $env_name" >&2; exit 1; }
}

ensure_env_by_module() {
  local module_key="$1"
  read -r env_name create_rel <<<"$(_env_spec "$module_key")"
  [[ "$env_name" == "UNKNOWN" ]] && { echo "[env] unknown module: $module_key" >&2; exit 1; }
  ensure_env "$env_name" "$ROOT/$create_rel"
}

run_in_env() {
  local env_name="$1"; shift
  micromamba run -n "$env_name" "$@"
}
