#!/usr/bin/env bash
set -euo pipefail
[[ "${DEBUG:-0}" == "1" ]] && set -x

# ---- Micromamba root autodetect ----
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

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
command -v micromamba >/dev/null || { echo "micromamba not found under $MAMBA_ROOT_PREFIX/bin"; exit 127; }
export GTDBTK_DATA_PATH="${GTDBTK_DATA_PATH:-$ROOT/refdata/gtdbtk/release226}"

FORCE_ENV_MODE="none"; FORCE_ENV_TARGET=""
MODE="full"; RESUME=1; FORCE=0; FROM=""; ONLY=""; SAMPLE=""
declare -a SETS=(); declare -a ARGS=()

usage(){ cat <<'USAGE'
Usage:
  bin/main.sh list
  bin/main.sh run <module|all> [--test|--resume|--force] [--from X] [--only X] [--sample ID] \
                               [--set k=v ...] [--force-env[=all|<module>|<env>]]
  bin/main.sh env list
  bin/main.sh env create <all|ingest|qc|assembly|binning|sra_tools|sim|plasmids>
  bin/main.sh env remove <ingest|qc|assembly|binning|sra_tools|sim>

Direct helpers:
  bin/main.sh plasmids_env
  bin/main.sh plasmids --sample ID [--cpus N] [--force|--resume]

Modules:
  ingest, qc, assembly, binning, sim, plasmids, all
USAGE
}

require_val(){ [[ $# -ge 2 && -n "${2:-}" && "${2:0:2}" != "--" ]] || { echo "error: $1 requires a value" >&2; exit 2; }; }

# shellcheck source=bin/lib/env.sh
source "$ROOT/bin/lib/env.sh"

merge(){ # $1 out_yaml, $2 module-defaults-yaml
  ensure_env_by_module "sra_tools"
  local setflags=(); for s in "${SETS[@]}"; do setflags+=(--set "$s"); done
  local m; m="$(mktemp)" || { echo "[err] failed to create temp file"; exit 2; }
  trap 'rm -f "$m"' RETURN
  cat "$2" > "$m"
  run_in_env "env_sra_tools" python "$ROOT/scripts/merge_params.py" \
    --module-defaults "$m" --defaults "$ROOT/config/defaults.yaml" --local "$ROOT/config/local.yaml" \
    "${setflags[@]}" --out "$1"
}

run_ingest_wrap(){ ensure_env_by_module "ingest"; local m p; m="$(mktemp)"; p="$(mktemp)"; trap 'rm -f "$m" "$p"' RETURN
  source "$ROOT/bin/modules/ingest.sh"; module_default_params > "$m"; merge "$p" "$m"
  export PARAMS_YAML="$p" ROOT MODE RESUME FORCE FROM ONLY SAMPLE
  run_in_env "env_ingest" bash -lc 'source "$ROOT/bin/modules/ingest.sh"; run_ingest "$@"' _ "${ARGS[@]:-}"; }

run_qc_wrap(){ ensure_env_by_module "qc"; local m p; m="$(mktemp)"; p="$(mktemp)"; trap 'rm -f "$m" "$p"' RETURN
  source "$ROOT/bin/modules/qc.sh"; module_default_params > "$m"; merge "$p" "$m"
  export PARAMS_YAML="$p" ROOT MODE RESUME FORCE FROM ONLY SAMPLE
  run_in_env "env_qc_core" bash -lc 'source "$ROOT/bin/modules/qc.sh"; run_qc "$@"' _ "${ARGS[@]:-}"; }

run_assembly_wrap(){ ensure_env_by_module "assembly"; local m p; m="$(mktemp)"; p="$(mktemp)"; trap 'rm -f "$m" "$p"' RETURN
  source "$ROOT/bin/modules/assembly.sh"; module_default_params > "$m"; merge "$p" "$m"
  export PARAMS_YAML="$p" ROOT MODE RESUME FORCE FROM ONLY SAMPLE
  run_in_env "env_assembly_core" bash -lc 'source "$ROOT/bin/modules/assembly.sh"; run_assembly "$@"' _ "${ARGS[@]:-}"; }

run_binning_wrap(){ ensure_env_by_module "binning"; local m p; m="$(mktemp)"; p="$(mktemp)"; trap 'rm -f "$m" "$p"' RETURN
  source "$ROOT/bin/modules/binning.sh"; module_default_params > "$m"; merge "$p" "$m"
  export PARAMS_YAML="$p" ROOT MODE RESUME FORCE FROM ONLY SAMPLE
  run_in_env "env_binning" bash -lc 'source "$ROOT/bin/modules/binning.sh"; run_binning "$@"' _ "${ARGS[@]:-}"; }

run_sim_wrap(){ export ROOT="$ROOT"; export CFG="${CFG:-$ROOT/config/sim.yaml}"
  [[ -s "$ROOT/bin/modules/sim.sh" ]] || { echo "[err] bin/modules/sim.sh not found"; exit 2; }
  bash "$ROOT/bin/modules/sim.sh"; }

run_plasmids_wrap(){
  [[ -s "$ROOT/bin/modules/plasmids.sh" ]] || { echo "[err] bin/modules/plasmids.sh not found"; exit 2; }
  export ROOT MODE RESUME FORCE FROM ONLY SAMPLE
  bash "$ROOT/bin/modules/plasmids.sh"
}

# ---- new: run all ----
run_all_wrap(){
  echo "[all] start MODE=$MODE RESUME=$RESUME FORCE=$FORCE SAMPLE=${SAMPLE:-}"
  run_ingest_wrap
  run_qc_wrap
  run_assembly_wrap
  run_binning_wrap
  # sim is optional and usually separate; keep out of default chain
  run_plasmids_wrap
  echo "[all] done"
}

env_list(){ printf "module\tenv\tcreator\n"; list_env_specs; }

env_create(){
  local target="${1:-}"; [[ -z "$target" ]] && { echo "[err] env create <target>"; exit 2; }
  if [[ "$target" == "all" ]]; then
    for m in ingest qc assembly binning sra_tools sim plasmids; do create_env_by_module "$m"; done
  else
    create_env_by_module "$target"
  fi
}

env_remove(){ local target="${1:-}"; [[ -z "$target" ]] && { echo "[err] env remove <target>"; exit 2; }; remove_env_by_module "$target"; }

cmd="${1:-}"; shift || true
case "$cmd" in
  list)
    ls -1 "$ROOT/bin/modules"/*.sh 2>/dev/null | xargs -n1 basename | sed 's/\.sh$//' || true
    ;;
  env)
    sub="${1:-}"; shift || true
    case "$sub" in
      list) env_list ;;
      create) env_create "${1:-}" ;;
      remove) env_remove "${1:-}" ;;
      *) usage; exit 2 ;;
    esac
    ;;
  plasmids_env)
    "$ROOT/env/create_env_plasmids.sh"; exit 0
    ;;
  plasmids)
    while [[ $# -gt 0 ]]; do
      case "$1" in
        --sample) require_val "$1" "${2:-}"; SAMPLE="$2"; shift ;;
        --sample=*) SAMPLE="${1#--sample=}";;
        --force) FORCE=1; RESUME=0 ;;
        --resume) RESUME=1; FORCE=0 ;;
        --cpus) require_val "$1" "${2:-}"; export CPUS="$2"; shift ;;
        --cpus=*) export CPUS="${1#--cpus=}";;
        *) ARGS+=("$1");;
      esac; shift || true
    done
    : "${SAMPLE:?need --sample or env SAMPLE}"
    run_plasmids_wrap
    ;;
  run)
    mod="${1:-}"; [[ -z "$mod" ]] && { usage; exit 2; }; shift || true
    while [[ $# -gt 0 ]]; do
      case "$1" in
        --test) MODE="test" ;;
        --resume) RESUME=1; FORCE=0 ;;
        --force) FORCE=1; RESUME=0 ;;
        --from)   require_val "$1" "${2:-}"; FROM="$2"; shift ;;
        --only)   require_val "$1" "${2:-}"; ONLY="$2"; shift ;;
        --sample) require_val "$1" "${2:-}"; SAMPLE="$2"; shift ;;
        --set)    require_val "$1" "${2:-}"; SETS+=("$2"); shift ;;
        --set=*)  SETS+=("${1#--set=}") ;;
        --fastq|--srr|--fasta) require_val "$1" "${2:-}"; ARGS+=("$1" "$2"); shift ;;
        --fastq=*|--srr=*|--fasta=*) ARGS+=("$1") ;;
        --force-env|--force-env=all) FORCE_ENV_MODE="all" ;;
        --force-env=*) FORCE_ENV_MODE="one"; FORCE_ENV_TARGET="${1#*=}" ;;
        *) ARGS+=("$1") ;;
      esac; shift || true
    done
    if [[ "$FORCE_ENV_MODE" == "one" ]]; then
      case "$FORCE_ENV_TARGET" in
        ingest)    FORCE_ENV_TARGET="env_ingest" ;;
        qc)        FORCE_ENV_TARGET="env_qc" ;;
        assembly)  FORCE_ENV_TARGET="env_assembly_core" ;;
        binning)   FORCE_ENV_TARGET="env_binning" ;;
        sra_tools) FORCE_ENV_TARGET="env_sra_tools" ;;
        sim)       FORCE_ENV_TARGET="sim_env" ;;
        plasmids)  FORCE_ENV_TARGET="env_plasmids" ;;
        *) : ;;
      esac
    fi
    export FORCE_ENV_MODE FORCE_ENV_TARGET
    case "$mod" in
      ingest)   run_ingest_wrap   ;;
      qc)       run_qc_wrap       ;;
      assembly) run_assembly_wrap ;;
      binning)  run_binning_wrap  ;;
      sim)      run_sim_wrap      ;;
      plasmids) run_plasmids_wrap ;;
      all)      : "${SAMPLE:?need --sample for plasmids in all}"; run_all_wrap ;;
      *) usage; exit 2 ;;
    esac
    ;;
  *) usage; exit 2 ;;
esac
