#!/usr/bin/env bash
set -euo pipefail
# env in: ROOT PARAMS_YAML SAMPLE MODE IN1 IN2 FORCE

OUTDIR="$ROOT/SRA/qc/fastp/$SAMPLE"
LOGS="$ROOT/logs"; EVENTS="$ROOT/api/status.jsonl"
mkdir -p "$OUTDIR" "$LOGS"

# cleanup on --force
if [[ "${FORCE:-0}" == "1" ]]; then
  rm -f "$OUTDIR/${SAMPLE}_trimmed_"*.fastq.gz "$OUTDIR/${SAMPLE}.fastp."* || true
fi

ts(){ date -u +"%Y-%m-%dT%H:%M:%SZ"; }
echo "{\"ts\":\"$(ts)\",\"module\":\"qc\",\"program\":\"fastp\",\"sample\":\"$SAMPLE\",\"phase\":\"start\"}" >> "$EVENTS"

# params
getp(){ micromamba run -n env_sra_tools python "$ROOT/scripts/merge_params.py" --defaults "$PARAMS_YAML" --get "$1"; }
LEN_REQ=$(getp modules.qc.fastp.length_required)
Q_PHRED=$(getp modules.qc.fastp.qualified_quality_phred)
UNQ_PCT=$(getp modules.qc.fastp.unqualified_percent_limit)
THREADS=$(getp global.threads)

# outputs
T1="$OUTDIR/${SAMPLE}_trimmed_1.fastq.gz"
T2="$OUTDIR/${SAMPLE}_trimmed_2.fastq.gz"
HTML="$OUTDIR/${SAMPLE}.fastp.html"
JSON="$OUTDIR/${SAMPLE}.fastp.json"

# log chosen inputs clearly
{
  echo "[fastp] MODE=$MODE SAMPLE=$SAMPLE"
  echo "[fastp] IN1=$IN1"
  echo "[fastp] IN2=${IN2:-}"
} >> "$LOGS/qc_fastp_${SAMPLE}.log"

# run fastp
if [[ -n "${IN2:-}" ]]; then
  micromamba run -n env_qc_fastp_fastqc_multiqc fastp \
    -i "$IN1" -I "$IN2" -o "$T1" -O "$T2" \
    -l "$LEN_REQ" -q "$Q_PHRED" -u "$UNQ_PCT" \
    -h "$HTML" -j "$JSON" -w "$THREADS" \
    >> "$LOGS/qc_fastp_${SAMPLE}.log" 2>&1
else
  micromamba run -n env_qc_fastp_fastqc_multiqc fastp \
    -i "$IN1" -o "$T1" \
    -l "$LEN_REQ" -q "$Q_PHRED" -u "$UNQ_PCT" \
    -h "$HTML" -j "$JSON" -w "$THREADS" \
    >> "$LOGS/qc_fastp_${SAMPLE}.log" 2>&1
  T2=""
fi

# event done
if [[ -n "$T2" ]]; then
  echo "{\"ts\":\"$(ts)\",\"module\":\"qc\",\"program\":\"fastp\",\"sample\":\"$SAMPLE\",\"phase\":\"done\",\"out\":{\"r1\":\"$T1\",\"r2\":\"$T2\"}}" >> "$EVENTS"
else
  echo "{\"ts\":\"$(ts)\",\"module\":\"qc\",\"program\":\"fastp\",\"sample\":\"$SAMPLE\",\"phase\":\"done\",\"out\":{\"r1\":\"$T1\"}}" >> "$EVENTS"
fi
