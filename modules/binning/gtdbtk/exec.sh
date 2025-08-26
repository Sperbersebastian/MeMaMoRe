#!/usr/bin/env bash
set -euo pipefail
BINS="$ROOT/SRA/binning/magscot/$SAMPLE/bins"
OUT="$ROOT/SRA/binning/gtdbtk/$SAMPLE"
mkdir -p "$OUT"

# if no refined bins, skip
shopt -s nullglob
files=("$BINS"/*.fa*); [[ ${#files[@]} -gt 0 ]] || { echo "[gtdbtk] no refined bins"; exit 0; }

# require ref data
: "${GTDBTK_DATA_PATH:?GTDBTK_DATA_PATH not set}"

micromamba run -n env_binning gtdbtk classify_wf \
  --genome_dir "$BINS" --out_dir "$OUT" --extension fa \
  --cpus 8 --skip_ani_screen || true

echo "[gtdbtk] -> $OUT"
