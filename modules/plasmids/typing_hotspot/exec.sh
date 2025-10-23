#!/usr/bin/env bash
# HOTSPOT host prediction (preprocess -> predict)
set -euo pipefail
ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

IN="$ROOT/SRA/plasmids/$SAMPLE/cluster/plasmids_derep.$SAMPLE.fasta"
OUT="$ROOT/SRA/plasmids/$SAMPLE/hotspot"; mkdir -p "$OUT"
[[ -s "$IN" ]] || { echo "[hotspot] no input: $IN"; exit 1; }

# Paths
HOTSPOT_DIR="${HOTSPOT_DIR:-$ROOT/external/HOTSPOT}"
HOTSPOT_DB="${HOTSPOT_DB:-$HOTSPOT_DIR/database}"
HOTSPOT_MODELS="${HOTSPOT_MODELS:-$HOTSPOT_DIR/models}"

[[ -d "$HOTSPOT_DIR" ]]   || { echo "[hotspot] repo missing: $HOTSPOT_DIR"; exit 2; }
[[ -d "$HOTSPOT_DB" ]]    || { echo "[hotspot] database missing: $HOTSPOT_DB"; exit 3; }
[[ -d "$HOTSPOT_MODELS" ]]|| { echo "[hotspot] models missing: $HOTSPOT_MODELS"; exit 4; }
command -v python >/dev/null || { echo "[hotspot] python not in PATH"; exit 5; }

DONE="$OUT/.done"
if [[ "$FORCE" == "1" ]]; then rm -f "$DONE"; fi
if [[ -s "$DONE" ]]; then echo "[hotspot] cached"; exit 0; fi

# Run HOTSPOT in repo working directory, then collect outputs into $OUT
pushd "$HOTSPOT_DIR" >/dev/null

# Step 1: preprocessing (encodes sequences)
python preprocessing.py \
  --fasta "$IN" \
  --database "$HOTSPOT_DB" \
  --model_path "$HOTSPOT_MODELS"

# Step 2: prediction
python hotspot.py

popd >/dev/null

# Collect outputs
if [[ -s "$HOTSPOT_DIR/results/host_lineage.tsv" ]]; then
  cp -f "$HOTSPOT_DIR/results/host_lineage.tsv" "$OUT/host_lineage.tsv"
  rsync -a "$HOTSPOT_DIR/results/" "$OUT/results/" >/dev/null 2>&1 || true
  date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
  echo "[hotspot] done -> $OUT/host_lineage.tsv"
else
  echo "[hotspot] prediction results not found in $HOTSPOT_DIR/results"; exit 6
fi
