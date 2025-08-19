#!/usr/bin/env bash
set -euo pipefail
OUTDIR="$ROOT/SRA/qc/fastqc/$SAMPLE"; LOGS="$ROOT/logs"; EVENTS="$ROOT/api/status.jsonl"
mkdir -p "$OUTDIR" "$LOGS"
ts(){ date -u +"%Y-%m-%dT%H:%M:%SZ"; }
echo "{\"ts\":\"$(ts)\",\"module\":\"qc\",\"program\":\"fastqc\",\"sample\":\"$SAMPLE\",\"phase\":\"start\"}" >> "$EVENTS"
if [[ -n "${IN2:-}" ]]; then IN="$IN1 $IN2"; else IN="$IN1"; fi
micromamba run -n env_qc_fastp_fastqc_multiqc fastqc -o "$OUTDIR" -t "$(micromamba run -n env_sra_tools python "$ROOT/scripts/merge_params.py" --defaults "$PARAMS_YAML" --get global.threads)" $IN \
  >"$LOGS/qc_fastqc_${SAMPLE}.log" 2>&1
echo "{\"ts\":\"$(ts)\",\"module\":\"qc\",\"program\":\"fastqc\",\"sample\":\"$SAMPLE\",\"phase\":\"done\"}" >> "$EVENTS"
