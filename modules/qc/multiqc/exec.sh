#!/usr/bin/env bash
set -euo pipefail
OUTDIR="$ROOT/SRA/qc/multiqc"; LOGS="$ROOT/logs"; EVENTS="$ROOT/api/status.jsonl"
mkdir -p "$OUTDIR" "$LOGS"
ts(){ date -u +"%Y-%m-%dT%H:%M:%SZ"; }
echo "{\"ts\":\"$(ts)\",\"module\":\"qc\",\"program\":\"multiqc\",\"phase\":\"start\"}" >> "$EVENTS"
micromamba run -n env_qc_fastp_fastqc_multiqc multiqc "$ROOT/SRA/qc" -o "$OUTDIR" \
  >"$LOGS/qc_multiqc.log" 2>&1
echo "{\"ts\":\"$(ts)\",\"module\":\"qc\",\"program\":\"multiqc\",\"phase\":\"done\",\"out\":{\"report\":\"$OUTDIR/multiqc_report.html\"}}" >> "$EVENTS"
