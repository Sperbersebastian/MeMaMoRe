#!/usr/bin/env bash
set -euo pipefail
[[ "${DEBUG:-0}" == "1" ]] && set -x

export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
command -v micromamba >/dev/null || { echo "micromamba not found under $MAMBA_ROOT_PREFIX/bin"; exit 127; }
export GTDBTK_DATA_PATH="${GTDBTK_DATA_PATH:-$ROOT/refdata/gtdbtk/release226}"

# force-env controls (default)
FORCE_ENV_MODE=none     # none|all|one
FORCE_ENV_TARGET=""

MODE="full"; RESUME=1; FORCE=0; FROM=""; ONLY=""; SAMPLE=""
declare -a SETS=()
declare -a ARGS=()

usage(){ cat <<EOF
Usage:
  bin/main.sh list
  bin/main.sh run <module> [--test|--resume|--force] [--from X] [--only X] [--sample ID] [--set k=v ...] [--force-env[=all|<module>|<env>]]
Modules:
  ingest, qc, assembly, binning
EOF
}

require_val(){ # $1=flag $2=next
  [[ $# -ge 2 && -n "${2:-}" && "${2:0:2}" != "--" ]] || { echo "error: $1 requires a value" >&2; exit 2; }
}

# Load env helper (after ROOT is set)
# shellcheck source=bin/lib/env.sh
source "$ROOT/bin/lib/env.sh"

merge(){ # $1 out, $2 module-defaults
  # Ensure the Python helper env exists
  ensure_env_by_module "sra_tools"

  # Build --set flags robustly
  local setflags=()
  for s in "${SETS[@]}"; do setflags+=(--set "$s"); done

  run_in_env "env_sra_tools" python "$ROOT/scripts/merge_params.py" \
    --module-defaults "$2" --defaults "$ROOT/config/defaults.yaml" \
    --local "$ROOT/config/local.yaml" "${setflags[@]}" --out "$1"
}

run_ingest_wrap(){
  source "$ROOT/bin/modules/ingest.sh"
  ensure_env_by_module "ingest"
  m="$(mktemp)"; module_default_params > "$m"
  p="$(mktemp)"; merge "$p" "$m"
  PARAMS_YAML="$p" ROOT="$ROOT" MODE="$MODE" RESUME="$RESUME" FORCE="$FORCE" FROM="$FROM" ONLY="$ONLY" SAMPLE="$SAMPLE" \
    run_in_env "env_ingest" bash -c 'run_ingest "$@"' _ run_ingest "${ARGS[@]:-}"
}

run_qc_wrap(){
  source "$ROOT/bin/modules/qc.sh"
  ensure_env_by_module "qc"
  m="$(mktemp)"; module_default_params > "$m"
  p="$(mktemp)"; merge "$p" "$m"
  PARAMS_YAML="$p" ROOT="$ROOT" MODE="$MODE" RESUME="$RESUME" FORCE="$FORCE" FROM="$FROM" ONLY="$ONLY" SAMPLE="$SAMPLE" \
    run_in_env "env_qc" bash -c 'run_qc "$@"' _ run_qc "${ARGS[@]:-}"
}

run_assembly_wrap(){
  source "$ROOT/bin/modules/assembly.sh"
  ensure_env_by_module "assembly"
  m="$(mktemp)"; module_default_params > "$m"
  p="$(mktemp)"; merge "$p" "$m"
  PARAMS_YAML="$p" ROOT="$ROOT" MODE="$MODE" RESUME="$RESUME" FORCE="$FORCE" FROM="$FROM" ONLY="$ONLY" SAMPLE="$SAMPLE" \
    run_in_env "env_assembly" bash -c 'run_assembly "$@"' _ run_assembly "${ARGS[@]:-}"
}

run_binning_wrap(){
  source "$ROOT/bin/modules/binning.sh"
  ensure_env_by_module "binning"
  m="$(mktemp)"; module_default_params > "$m"
  p="$(mktemp)"; merge "$p" "$m"
  export PARAMS_YAML="$p" ROOT="$ROOT" MODE="$MODE" RESUME="$RESUME" FORCE="$FORCE" FROM="$FROM" ONLY="$ONLY" SAMPLE="$SAMPLE"
  run_in_env "env_binning" bash -lc 'source "$ROOT/bin/modules/binning.sh"; run_binning "$@"' _ "${ARGS[@]:-}"
}


cmd="${1:-}"; shift || true

case "$cmd" in
  list)
    ls -1 "$ROOT/bin/modules"/*.sh 2>/dev/null | xargs -n1 basename | sed 's/\.sh$//' || true
    ;;
  run)
    mod="${1:-}"; [[ -z "$mod" ]] && { usage; exit 2; }; shift || true
    declare -a ARGS=()
    while [[ $# -gt 0 ]]; do
      case "$1" in
        --test) MODE="test";;
        --resume) RESUME=1; FORCE=0;;
        --force) FORCE=1; RESUME=0;;
        --from)   require_val "$1" "${2:-}"; FROM="$2"; shift;;
        --only)   require_val "$1" "${2:-}"; ONLY="$2"; shift;;
        --sample) require_val "$1" "${2:-}"; SAMPLE="$2"; shift;;
        --set)    require_val "$1" "${2:-}"; SETS+=("$2"); shift;;
        --set=*)  SETS+=("${1#--set=}");;

        # IO shortcuts passed to module
        --fastq|--srr|--fasta) require_val "$1" "${2:-}"; ARGS+=("$1" "$2"); shift;;
        --fastq=*|--srr=*|--fasta=*) ARGS+=("$1");;

        # env forcing
        --force-env)            FORCE_ENV_MODE=all ;;
        --force-env=all)        FORCE_ENV_MODE=all ;;
        --force-env=*)
          FORCE_ENV_MODE=one
          FORCE_ENV_TARGET="${1#*=}"    # module key or env name
          ;;
        *) ARGS+=("$1");;
      esac
      shift || true
    done

    # Translate module key -> env name for single-target force
    if [[ "$FORCE_ENV_MODE" == "one" ]]; then
      case "$FORCE_ENV_TARGET" in
        ingest)   FORCE_ENV_TARGET="env_ingest" ;;
        qc)       FORCE_ENV_TARGET="env_qc" ;;
        assembly) FORCE_ENV_TARGET="env_assembly" ;;
        binning)  FORCE_ENV_TARGET="env_binning" ;;
        sra_tools)FORCE_ENV_TARGET="env_sra_tools" ;;
      esac
    fi

    # Export for env helper
    export FORCE_ENV_MODE FORCE_ENV_TARGET

    case "$mod" in
      ingest)   run_ingest_wrap   ;;
      qc)       run_qc_wrap       ;;
      assembly) run_assembly_wrap ;;
      binning)  run_binning_wrap  ;;
      *) usage; exit 2;;
    esac
    ;;
  *) usage; exit 2;;
esac
