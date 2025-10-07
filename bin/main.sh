#!/usr/bin/env bash
set -euo pipefail

export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

MODE="full"; RESUME=1; FORCE=0; SETS=(); FROM=""; ONLY=""; SAMPLE=""

usage(){ cat <<EOF
Usage:
  bin/main.sh list
  bin/main.sh run ingest --fastq <manifest> [--set k=v ...]
  bin/main.sh run qc [--test] [--sample ID] [--set k=v ...] [--force]
EOF
}

merge(){  # $1 out, $2 module-defaults
  micromamba run -n env_sra_tools python "$ROOT/scripts/merge_params.py" \
    --module-defaults "$2" --defaults "$ROOT/config/defaults.yaml" \
    --local "$ROOT/config/local.yaml" "${SETS[@]/#/--set }" --out "$1"
}

run_ingest_wrap(){
  # shellcheck source=/dev/null
  source "$ROOT/bin/modules/ingest.sh"
  m=/tmp/ingest_defaults.yaml; module_default_params > "$m"
  p=/tmp/ingest_params.yaml; merge "$p" "$m"
  PARAMS_YAML="$p" ROOT="$ROOT" MODE="$MODE" RESUME="$RESUME" FORCE="$FORCE" FROM="$FROM" ONLY="$ONLY" SAMPLE="$SAMPLE" run_ingest "$@"
}

run_qc_wrap(){
  # shellcheck source=/dev/null
  source "$ROOT/bin/modules/qc.sh"
  m=/tmp/qc_defaults.yaml; module_default_params > "$m"
  p=/tmp/qc_params.yaml; merge "$p" "$m"
  PARAMS_YAML="$p" ROOT="$ROOT" MODE="$MODE" RESUME="$RESUME" FORCE="$FORCE" FROM="$FROM" ONLY="$ONLY" SAMPLE="$SAMPLE" run_qc "$@"
}

run_assembly_wrap(){
  # shellcheck source=/dev/null
  source "$ROOT/bin/modules/assembly.sh"
  m=/tmp/assembly_defaults.yaml; module_default_params > "$m"
  p=/tmp/assembly_params.yaml; merge "$p" "$m"
  PARAMS_YAML="$p" ROOT="$ROOT" MODE="$MODE" RESUME="$RESUME" FORCE="$FORCE" FROM="$FROM" ONLY="$ONLY" SAMPLE="$SAMPLE" run_assembly "$@"
}

cmd="${1:-}"; shift || true

case "$cmd" in
  list)
    ls -1 "$ROOT/bin/modules" | sed 's/\.sh$//' || true
    ;;
  run)
    # 1) get module
    mod="${1:-}"; [[ -z "$mod" ]] && { usage; exit 2; }; shift || true
    # 2) parse flags for this run
    while [[ $# -gt 0 ]]; do
      case "$1" in
        --test) MODE="test";;
        --resume) RESUME=1; FORCE=0;;
        --force) FORCE=1; RESUME=0;;
        --from) FROM="${2:-}"; shift;;
        --only) ONLY="${2:-}"; shift;;
        --sample) SAMPLE="${2:-}"; shift;;
        --set) SETS+=("${2:-}"); shift;;
        --set=*) SETS+=("${1#--set=}");;
        --fastq|--srr|--fasta) ARGS+=("$1" "$2"); shift;;
        --fastq=*|--srr=*|--fasta=*) ARGS+=("$1");;
        *) ARGS+=("$1");;
      esac
      shift || true
    done
    # 3) dispatch
    case "$mod" in
      ingest)   run_ingest_wrap   "${ARGS[@]:-}";;
      qc)       run_qc_wrap       "${ARGS[@]:-}";;
      assembly) run_assembly_wrap "${ARGS[@]:-}";;
      *) usage; exit 2;;
    esac
    ;;
  *) usage; exit 2;;
esac
