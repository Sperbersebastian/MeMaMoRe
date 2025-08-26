#!/usr/bin/env bash
set -euo pipefail
# env: ROOT SAMPLE
BINS="$ROOT/SRA/binning/magscot/$SAMPLE/bins"
RR="$ROOT/SRA/binning/barrnap/$SAMPLE"
TR="$ROOT/SRA/binning/trnascan/$SAMPLE"
mkdir -p "$RR" "$TR"
[[ -d "$BINS" ]] || { echo "[rRNA/tRNA] no refined bins"; exit 0; }

for f in "$BINS"/*.fa*; do
  base=$(basename "$f"); base="${base%.*}"
  micromamba run -n env_binning barrnap --threads 8 < "$f" > "$RR/${base}.gff" 2>/dev/null || true
  micromamba run -n env_binning tRNAscan-SE -o "$TR/${base}.tsv" -q "$f" 2>/dev/null || true
done
echo "[rRNA/tRNA] -> $RR , $TR"
