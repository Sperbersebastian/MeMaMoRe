#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"

MODE="full"
RESUME=1
FORCE=0
SETS=()
FROM=""
ONLY=""
SAMPLE=""

usage() {
  cat <<EOF
Usage:
  bin/main.sh prepare-envs
  bin/main.sh list
  bin/main.sh run <module> [--test] [--set k=v ...] [--from X] [--only X] [--sample ID] [--resume|--force]
  bin/main.sh all [--test] [--set k=v ...] [--from X] [--sample ID] [--resume|--force]
EOF
}

# write merged YAML to $1; module defaults at $2
merge() {
  local out="$1" mdefs="$2"
  micromamba run -n env_sra_tools python "$ROOT/scripts/merge_params.py" \
    --module-defaults "$mdefs" \
    --defaults "$ROOT/config/defaults.yaml" \
    --local "$ROOT/config/local.yaml" \
    "${SETS[@]/#/--set }" \
    --out "$out"
}

run_ingest() {
  # shellcheck source=/dev/null
  source "$ROOT/bin/modules/ingest.sh"
  local m="/tmp/ingest_defaults.yaml"
  module_default_params > "$m"
  local p="/tmp/ingest_params.yaml"
  merge "$p" "$m"
  PARAMS_YAML="$p" ROOT="$ROOT" MODE="$MODE" RESUME="$RESUME" FORCE="$FORCE" FROM="$FROM" ONLY="$ONLY" SAMPLE="$SAMPLE" \
    run_ingest "$@"
}

list() {
  echo "Modules:"
  ls -1 "$ROOT/bin/modules" | sed 's/\.sh$//' || true
}

prepare_envs() {
  echo "Create micromamba envs with scripts in envs/. (stub)"
}

# ---- parse global flags after subcommand ----
cmd="${1:-}"; shift || true
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
    *) break;;
  esac
  shift || true
done

case "$cmd" in
  prepare-envs) prepare_envs;;
  list) list;;
  run)
    mod="${1:-}"; shift || true
    [[ -z "$mod" ]] && { usage; exit 2; }
    case "$mod" in
      ingest) run_ingest "$@";;
      *) echo "Unknown module: $mod"; list; exit 2;;
    esac
    ;;
  all)
    run_ingest "$@"
    ;;
  *) usage; exit 2;;
esac
