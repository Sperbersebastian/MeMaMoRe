#!/usr/bin/env bash
set -euo pipefail

: "${ROOT:?}"; : "${SAMPLE:?}"; : "${PARAMS_YAML:?}"
: "${IN1:?}"    # R1 path
IN2="${IN2:-}"  # optional R2 path

OUTDIR="$ROOT/SRA/qc/fastqc/$SAMPLE"
LOGS="$ROOT/logs"
EVENTS="$ROOT/api/status.jsonl"
mkdir -p "$OUTDIR" "$LOGS" "$(dirname "$EVENTS")"

ts(){ date -u +"%Y-%m-%dT%H:%M:%SZ"; }

echo "{\"ts\":\"$(ts)\",\"module\":\"qc\",\"program\":\"fastqc\",\"sample\":\"$SAMPLE\",\"phase\":\"start\"}" >> "$EVENTS"

# get threads once
THREADS="$(micromamba run -n env_sra_tools python "$ROOT/scripts/merge_params.py" --defaults "$PARAMS_YAML" --get global.threads)"

# build input array safely
inputs=("$IN1")
[[ -n "$IN2" ]] && inputs+=("$IN2")

# run fastqc in the correct env
micromamba run -n env_qc_core fastqc \
  -o "$OUTDIR" -t "$THREADS" \
  "${inputs[@]}" \
  >"$LOGS/qc_fastqc_${SAMPLE}.log" 2>&1

echo "{\"ts\":\"$(ts)\",\"module\":\"qc\",\"program\":\"fastqc\",\"sample\":\"$SAMPLE\",\"phase\":\"done\"}" >> "$EVENTS"
