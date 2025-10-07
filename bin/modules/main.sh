#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MODE="full"; RESUME=1; FORCE=0
SETS=(); FROM=""; ONLY=""; SAMPLE=""

usage(){ cat <<EOF
Usage:
  bin/main.sh prepare-envs
  bin/main.sh list
  bin/main.sh run <module> [--test] [--set k=v ...] [--from X] [--only X] [--sample ID] [--resume|--force]
  bin/main.sh all [--test] [--set k=v ...] [--from X] [--sample ID] [--resume|--force]
EOF
}

merge_params(){  # writes merged YAML to $1
  local out="$1" mdefs="$2"
  python "$ROOT/scripts/merge_params.py" \
    --module-defaults "$mdefs" \
    --defaults "$ROOT/config/defaults.yaml" \
    --local "$ROOT/config/local.yaml" \
    "${SETS[@]/#/--set }" \
    --out "$out"
}

timestamp(){ date -u +"%Y%m%dT%H%M%SZ"; }

prepare_envs(){
  echo "(stub) create env scripts here (envs/create_*.sh) and run them."
}

list(){
  echo "Modules:"; ls -1 "$ROOT/bin/modules" | sed 's/\.sh$//' || true
}

run_module(){  # $1=module
  local mod="$1"; shift || true
  local mfile="$ROOT/bin/modules/${mod}.sh"
  [[ -f "$mfile" ]] || { echo "Unknown module: $mod"; exit 2; }
  source "$mfile"

  local ts="$(timestamp)" rundir="$ROOT/runs/$ts"; mkdir -p "$rundir"
  local mdefs="/tmp/${mod}_defaults.yaml"
  module_default_params > "$mdefs"

  local params="$rundir/params.${mod}.yaml"
  merge_params "$params" "$mdefs"

  # Write plan skeleton for GUI
  echo "{\"ts\":\"$ts\",\"module\":\"$mod\",\"mode\":\"$MODE\",\"sample\":\"${SAMPLE:-all}\"}" > "$rundir/plan.json"
  echo "{\"ts\":\"$ts\",\"cmd\":\"run\",\"module\":\"$mod\"}" > "$rundir/run.json"

  PARAMS_YAML="$params" MODE="$MODE" RESUME="$RESUME" FORCE="$FORCE" FROM="$FROM" ONLY="$ONLY" SAMPLE="$SAMPLE" ROOT="$ROOT" \
    "run_${mod}"
}

cmd="${1:-}"; shift || true
while [[ $# -gt 0 ]]; do
  case "$1" in
    --test) MODE="test";;
    --resume) RESUME=1; FORCE=0;;
    --force) FORCE=1; RESUME=0;;
    --from) FROM="$2"; shift;;
    --only) ONLY="$2"; shift;;
    --sample) SAMPLE="$2"; shift;;
    --set) SETS+=("$2"); shift;;
    --set=*) SETS+=("${1#--set=}");;
    *) break;;
  esac; shift || true
done

case "$cmd" in
  prepare-envs) prepare_envs;;
  list) list;;
  run) mod="${1:-}"; [[ -z "$mod" ]] && usage || run_module "$mod";;
  all) for m in ingest qc assembly; do run_module "$m"; done;;
  *) usage; exit 2;;
esac
