#!/usr/bin/env bash
set -euo pipefail
# env: ROOT SAMPLE

BINS="$ROOT/SRA/binning/magscot/$SAMPLE/bins"
RR="$ROOT/SRA/binning/barrnap/$SAMPLE"
TR="$ROOT/SRA/binning/trnascan/$SAMPLE"
LOG="$ROOT/logs/binning_qc_rrna_trna_${SAMPLE}.log"

# reset outputs
rm -rf "$RR" "$TR" "$LOG"
mkdir -p "$RR" "$TR" "$(dirname "$LOG")"

[[ -d "$BINS" ]] || { echo "[rRNA/tRNA] no refined bins" | tee -a "$LOG"; exit 0; }

for f in "$BINS"/*.fa*; do
  base=$(basename "$f"); base="${base%.*}"
  # reset per-file outputs
  rm -f "$RR/${base}.gff" "$TR/${base}.tsv"
  {
    echo "[barrnap] $base"
    micromamba run -n env_binning barrnap --threads 8 < "$f" > "$RR/${base}.gff"
    echo "[tRNAscan-SE] $base"
    micromamba run -n env_binning tRNAscan-SE -B -o "$TR/${base}.tsv" -q "$f"
  } >>"$LOG" 2>&1 || true
done

echo "[rRNA/tRNA] -> $RR , $TR" | tee -a "$LOG"
