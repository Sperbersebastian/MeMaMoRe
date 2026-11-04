#!/usr/bin/env bash
# modules/args/deeparg/exec.sh
# Run DeepARG (LS model) once on the central gene calls (genes.ffn)
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

OUT="$ROOT/SRA/args/$SAMPLE/deeparg"; mkdir -p "$OUT"
DONE="$OUT/.done"
[[ "$FORCE" == "1" ]] && rm -f "$DONE"
[[ -s "$DONE" ]] && { echo "[deeparg] cached"; exit 0; }

command -v deeparg >/dev/null 2>&1 || { echo "[deeparg] deeparg not in PATH"; exit 3; }

# Auto-download DeepARG model/database on first run
DEEPARG_DATA="${DEEPARG_DATA:-$ROOT/refdata/deeparg}"
export DEEPARG_DATA
# DeepARG ships model/v2/{model,metadata}_{LS,SS}.pkl (there is no plain metadata.pkl;
# checking for it re-downloaded the 1.9 GB data on every run)
if [[ ! -s "$DEEPARG_DATA/model/v2/metadata_LS.pkl" || ! -s "$DEEPARG_DATA/model/v2/model_LS.pkl" ]]; then
  echo "[deeparg] downloading model & database -> $DEEPARG_DATA"
  mkdir -p "$DEEPARG_DATA"
  deeparg download_data -o "$DEEPARG_DATA" || { echo "[deeparg] DB download failed"; exit 5; }
fi

# ---- input: central gene calls (DeepARG-LS is trained on full-length genes) ----
GENES="$ROOT/SRA/assemblies/genes/$SAMPLE/genes.ffn"
[[ -s "$GENES" ]] || { echo "[deeparg] missing $GENES (run: bin/main.sh run assembly --sample $SAMPLE)"; exit 2; }

# Every gene is scanned once; contig/MAG/plasmid/virus context is assigned
# afterwards from the gene's contig ID (args summary / arg_context.tsv).
CTX_DIR="$OUT/genes"; mkdir -p "$CTX_DIR"
OUTPFX="$CTX_DIR/genes"
if [[ -s "${OUTPFX}.mapping.ARG" && "$FORCE" != "1" && "${OUTPFX}.mapping.ARG" -nt "$GENES" ]]; then
  echo "[deeparg] cached -> ${OUTPFX}.mapping.ARG"
else
  echo "[deeparg] scanning $(grep -c '^>' "$GENES") genes: $GENES"
  deeparg predict \
    --model LS \
    --type nucl \
    -i "$GENES" \
    -o "$OUTPFX" \
    -d "$DEEPARG_DATA" \
    --arg-alignment-identity 50 \
    --arg-alignment-evalue 1e-10 \
    --arg-num-alignments-per-entry 1000 \
    2>&1
  [[ -f "${OUTPFX}.mapping.ARG" ]] || { echo "[deeparg] ERROR: no ${OUTPFX}.mapping.ARG"; exit 4; }
fi
echo "[deeparg] $(tail -n +2 "${OUTPFX}.mapping.ARG" | wc -l) ARG predictions"

date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
echo "[deeparg] done -> $OUT"
