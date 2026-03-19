#!/usr/bin/env bash
# Env-Manager + Micromamba-Init
set -euo pipefail

# ---------- Repo-Root robust bestimmen ----------
ROOT="${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../" && pwd)}"
have_layout(){ [[ -d "$1/bin" && -d "$1/envs" ]]; }

if ! have_layout "$ROOT"; then
  try="$(cd "$ROOT/.." && pwd)"
  have_layout "$try" && ROOT="$try"
fi

if ! have_layout "$ROOT"; then
  if command -v git >/dev/null 2>&1; then
    top="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel 2>/dev/null || true)"
    [[ -n "$top" ]] && have_layout "$top" && ROOT="$top"
  fi
fi
export ROOT

# ---------- Micromamba ----------
if [[ -z "${MAMBA_ROOT_PREFIX:-}" ]]; then
  if [[ -d "$HOME/micromamba/envs" ]]; then
    export MAMBA_ROOT_PREFIX="$HOME/micromamba"
  elif [[ -d "$HOME/.local/share/mamba/envs" ]]; then
    export MAMBA_ROOT_PREFIX="$HOME/.local/share/mamba"
  else
    export MAMBA_ROOT_PREFIX="$HOME/micromamba"
  fi
fi
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
command -v micromamba >/dev/null 2>&1 || { echo "[env] micromamba not found at $MAMBA_ROOT_PREFIX/bin"; exit 127; }
eval "$(micromamba shell hook --shell=bash)"

# ---------- Force-Flags (optional) ----------
FORCE_ENV_MODE="${FORCE_ENV_MODE:-none}"     # none|all|one
FORCE_ENV_TARGET="${FORCE_ENV_TARGET:-}"     # env- oder modul-name

# =========================
# Helpers
# =========================
_env_name(){
  case "$1" in
    ingest)                 echo "env_ingest" ;;
    qc)                     echo "env_qc" ;;
    assembly)               echo "env_assembly" ;;
    binning)                echo "env_binning" ;;
    sra_tools)              echo "env_sra_tools" ;;
    sim)                    echo "sim_env" ;;
    mapping|mapping_coverm) echo "env_mapping_coverm" ;;
    coverm)                 echo "env_mapping_coverm" ;;
    plasmids)               echo "env_plasmids" ;;
    # --- Plasmid-Submodule ---
    metaplasmidspades)      echo "env_assembly_core" ;;
    viralverify)            echo "viralverify_env" ;;
    plasme)                 echo "plasme_env" ;;
    genomad)                echo "genomad_env" ;;
    mobrecon)               echo "mobsuite_env" ;;
    union_cluster)          echo "plasmids_core" ;;
    typing_hotspot)         echo "hotspot_env" ;;
    typer)                  echo "mobsuite_env" ;;
    cluster)                echo "stampede_env" ;;
    # --- Virus-Submodule ---
    viruses)                echo "env_viruses" ;;
    metaviralspades)        echo "env_assembly_core" ;;
    virsorter2)             echo "virsorter2_env" ;;
    vibrant)                echo "vibrant_env" ;;
    union_virus)            echo "plasmids_core" ;;
    checkv)                 echo "checkv_env" ;;
    # --- ARG-Submodule ---
    args)                   echo "env_args" ;;
    deeparg)                echo "deeparg_env" ;;
    rgi)                    echo "rgi_env" ;;
    amrplusplus)            echo "amrplusplus_env" ;;
    summary_args)           echo "deeparg_env" ;;
    *)                      echo "UNKNOWN" ;;
  esac
}

# Creator nur unter envs/ suchen (kein env/)
_resolve_creator(){
  local mod="$1"
  local exact="$ROOT/envs/create_env_${mod}.sh"
  if [[ -s "$exact" ]]; then echo "$exact"; return 0; fi

  local first
  first="$(ls -1 "$ROOT"/envs/create_env_"$mod"*.sh 2>/dev/null | head -n1 || true)"
  [[ -n "$first" ]] && { echo "$first"; return 0; }

  case "$mod" in
    metaplasmidspades) echo "$ROOT/envs/create_env_assembly.sh"; return 0 ;;
    # All plasmid substeps are provisioned by the umbrella creator
    viralverify|plasme|genomad|mobrecon|union_cluster|typing_hotspot|coverm|cluster|typer)
      [[ -s "$ROOT/envs/create_env_plasmids.sh" ]] && { echo "$ROOT/envs/create_env_plasmids.sh"; return 0; }
      ;;
    # Virus substeps: metaviralspades reuses assembly env, others have own creators or reuse plasmids
    metaviralspades) echo "$ROOT/envs/create_env_assembly.sh"; return 0 ;;
    union_virus)
      [[ -s "$ROOT/envs/create_env_plasmids.sh" ]] && { echo "$ROOT/envs/create_env_plasmids.sh"; return 0; }
      ;;
    # ARG substeps: summary_args reuses deeparg env
    summary_args)
      [[ -s "$ROOT/envs/create_env_deeparg.sh" ]] && { echo "$ROOT/envs/create_env_deeparg.sh"; return 0; }
      ;;
  esac
  return 1
}

_list_env_names(){ micromamba env list | awk 'NR>2{print $1}' | sed '/^$/d'; }

_env_exists(){
  local name="$1"
  if   [[ "$name" == "env_qc" ]]; then
    _list_env_names | grep -qx env_qc_core && _list_env_names | grep -qx env_qc_mqc
  elif [[ "$name" == "env_assembly" ]]; then
    _list_env_names | grep -qx env_assembly_core && _list_env_names | grep -qx env_mapping_coverm
  elif [[ "$name" == "env_plasmids" ]]; then
    _list_env_names | grep -qx env_assembly_core &&
    _list_env_names | grep -qx plasmids_core   &&
    _list_env_names | grep -qx genomad_env     &&
    _list_env_names | grep -qx plasme_env      &&
    _list_env_names | grep -qx viralverify_env &&
    _list_env_names | grep -qx mobsuite_env    &&
    _list_env_names | grep -qx hotspot_env     &&
    _list_env_names | grep -qx stampede_env
  elif [[ "$name" == "env_viruses" ]]; then
    _list_env_names | grep -qx env_assembly_core &&
    _list_env_names | grep -qx viralverify_env   &&
    _list_env_names | grep -qx genomad_env       &&
    _list_env_names | grep -qx virsorter2_env    &&
    _list_env_names | grep -qx vibrant_env       &&
    _list_env_names | grep -qx checkv_env
  elif [[ "$name" == "env_args" ]]; then
    _list_env_names | grep -qx deeparg_env      &&
    _list_env_names | grep -qx rgi_env           &&
    _list_env_names | grep -qx amrplusplus_env
  else
    _list_env_names | grep -qx "$name"
  fi
}

_recreate_env(){
  local env_name="$1" create_script="$2"
  echo "[env] recreate $env_name via $create_script"
  micromamba env remove -n "$env_name" -y || true
  bash "$create_script"
}

# =========================
# Public API
# =========================
list_env_specs(){
  for m in ingest qc assembly binning sra_tools sim mapping_coverm coverm plasmids \
           viralverify plasme genomad mobrecon typer cluster union_cluster typing_hotspot \
           virsorter2 vibrant checkv \
           deeparg rgi amrplusplus summary_args; do
    local env create
    env="$(_env_name "$m")"
    if create="$(_resolve_creator "$m")"; then :; else
      create="(missing: $ROOT/envs/create_env_${m}.sh)"
    fi
    printf "%s\t%s\t%s\n" "$m" "$env" "$create"
  done
}

ensure_env(){
  local env_name="$1" create_script="$2"
  local force="no"
  [[ "$FORCE_ENV_MODE" == "all" ]] && force="yes"
  if [[ "$FORCE_ENV_MODE" == "one" && -n "${FORCE_ENV_TARGET:-}" ]]; then
    if [[ "$FORCE_ENV_TARGET" == "$env_name" || "$(_env_name "$FORCE_ENV_TARGET")" == "$env_name" ]]; then
      force="yes"
    fi
  fi

  case "$env_name" in
    env_qc|env_assembly|env_plasmids)
      if [[ "$force" == "yes" ]]; then
        echo "[env] recreate $env_name via $create_script"; bash "$create_script"
      else
        if ! _env_exists "$env_name"; then
          echo "[env] create $env_name via $create_script"; bash "$create_script"
        fi
      fi
      _env_exists "$env_name" || { echo "[env] failed to create $env_name"; exit 1; }
      ;;
    *)
      if [[ "$force" == "yes" ]]; then
        _recreate_env "$env_name" "$create_script"
      elif ! _env_exists "$env_name"; then
        echo "[env] create $env_name via $create_script"; bash "$create_script"
      fi
      _env_exists "$env_name" || { echo "[env] failed to create $env_name" >&2; exit 1; }
      ;;
  esac
}

ensure_env_by_module(){
  local mod="$1"
  local env create
  env="$(_env_name "$mod")"
  [[ "$env" == "UNKNOWN" ]] && { echo "[env] unknown module: $mod" >&2; exit 1; }
  if ! create="$(_resolve_creator "$mod")"; then
    echo "[env] creator script not found for module '$mod'. Expected:"
    echo "  $ROOT/envs/create_env_${mod}.sh"
    exit 1
  fi
  ensure_env "$env" "$create"
}

create_env_by_module(){
  local mod="$1"
  local env create
  env="$(_env_name "$mod")"
  [[ "$env" == "UNKNOWN" ]] && { echo "[env] unknown module: $mod" >&2; exit 1; }
  if ! create="$(_resolve_creator "$mod")"; then
    echo "[env] creator script not found for module '$mod'. Expected:"
    echo "  $ROOT/envs/create_env_${mod}.sh"
    exit 1
  fi

  case "$env" in
    env_qc|env_assembly|env_plasmids)
      echo "[env] recreate $env via $create"
      bash "$create"
      _env_exists "$env" || { echo "[env] failed to create $env"; exit 1; }
      ;;
    *)
      _recreate_env "$env" "$create"
      ;;
  esac
}

env_name_for_module(){
  local mod="$1"
  local env; env="$(_env_name "$mod")"
  [[ "$env" == "UNKNOWN" ]] && { echo "[env] unknown module: $mod" >&2; return 1; }
  echo "$env"
}

remove_env_by_module(){
  local mod="$1"
  local env; env="$(_env_name "$mod")"
  [[ "$env" == "UNKNOWN" ]] && { echo "[env] unknown module: $mod" >&2; exit 1; }

  case "$env" in
    env_qc)
      micromamba env remove -n env_qc_core -y || true
      micromamba env remove -n env_qc_mqc  -y || true
      ;;
    env_assembly)
      micromamba env remove -n env_assembly_core -y || true
      micromamba env remove -n env_mapping_coverm -y || true
      ;;
    env_plasmids)
      for e in env_assembly_core plasmids_core genomad_env plasme_env viralverify_env mobsuite_env hotspot_env stampede_env; do
        micromamba env remove -n "$e" -y || true
      done
      ;;
    *)
      micromamba env remove -n "$env" -y || true
      ;;
  esac
}

run_in_env(){ local env_name="$1"; shift; micromamba run -n "$env_name" "$@"; }

require_cmd(){ for c in "$@"; do command -v "$c" >/dev/null 2>&1 || { echo "[dep] missing: $c"; return 1; }; done; }
