#!/usr/bin/env bash
set -euo pipefail
# env: ROOT SAMPLE CONTIGS BAM [BAMS]
# BAMS (space-separated, binning group >= 2 samples): one depth column per BAM
# (multi-sample coverage for MetaBAT2/Binny); otherwise the sample's own BAM.
OUT="$ROOT/SRA/binning/coverage/$SAMPLE"
LOG="$ROOT/logs/binning_coverage_${SAMPLE}.log"
mkdir -p "$OUT"
DEP="$OUT/${SAMPLE}.depth.txt"

# Use the BAM as-is; sort/index defensively
SORTED_BAM="$OUT/${SAMPLE}.sorted.bam"
if [[ -s "$BAM" ]]; then
  if [[ -s "$BAM.bai" ]]; then
    cp -f "$BAM" "$SORTED_BAM" && cp -f "$BAM.bai" "$SORTED_BAM.bai"
  else
    micromamba run -n env_binning samtools sort -@ 4 "$BAM" -o "$SORTED_BAM"
    micromamba run -n env_binning samtools index "$SORTED_BAM"
  fi
else
  echo "[coverage] ERROR: BAM missing for $SAMPLE ($BAM)" | tee -a "$LOG"; exit 2
fi

# JGI without referenceFasta to avoid header mismatches
read -r -a DEPTH_BAMS <<<"${BAMS:-$SORTED_BAM}"
micromamba run -n env_binning jgi_summarize_bam_contig_depths \
  --outputDepth "$DEP" \
  --minMapQual 10 \
  --percentIdentity 97 \
  "${DEPTH_BAMS[@]}" >"$LOG" 2>&1
echo "[coverage] ${#DEPTH_BAMS[@]} sample(s) in depth table"

[[ -s "$DEP" ]] || { echo "[coverage] ERROR: depth empty for $SAMPLE" | tee -a "$LOG"; exit 2; }
echo "[coverage] -> $DEP"
