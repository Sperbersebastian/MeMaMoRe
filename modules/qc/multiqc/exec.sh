#!/usr/bin/env bash
set -euo pipefail

: "${ROOT:?}"; : "${SAMPLE:?}"

INROOT="$ROOT/SRA/qc"                  # where reports live
OUTDIR="$INROOT/multiqc/$SAMPLE"       # output dir for MultiQC
LOGS="$ROOT/logs"
EVENTS="$ROOT/api/status.jsonl"
mkdir -p "$OUTDIR" "$LOGS" "$(dirname "$EVENTS")"

ts(){ date -u +"%Y-%m-%dT%H:%M:%SZ"; }
echo "{\"ts\":\"$(ts)\",\"module\":\"qc\",\"program\":\"multiqc\",\"sample\":\"$SAMPLE\",\"phase\":\"start\"}" >> "$EVENTS"

# Point MultiQC to the sample’s QC subfolders (don’t exec a directory!)
micromamba run -n env_qc_mqc multiqc \
  -o "$OUTDIR" \
  "$INROOT/fastqc/$SAMPLE" \
  "$INROOT/fastp/$SAMPLE" \
  >"$LOGS/qc_multiqc_${SAMPLE}.log" 2>&1

echo "{\"ts\":\"$(ts)\",\"module\":\"qc\",\"program\":\"multiqc\",\"sample\":\"$SAMPLE\",\"phase\":\"done\"}" >> "$EVENTS"
