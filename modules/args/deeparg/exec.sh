#!/usr/bin/env bash
# modules/args/deeparg/exec.sh
# Run DeepARG on multiple genomic contexts: contigs, MAGs, plasmids, viruses
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
if [[ ! -d "$DEEPARG_DATA" || ! -f "$DEEPARG_DATA/model/v2/metadata.pkl" ]]; then
  echo "[deeparg] downloading model & database -> $DEEPARG_DATA"
  mkdir -p "$DEEPARG_DATA"
  deeparg download_data -o "$DEEPARG_DATA" || { echo "[deeparg] DB download failed"; exit 5; }
fi

# ---- helper: run DeepARG on one FASTA ----
_run_deeparg() {
  local fasta="$1" context="$2" label="$3"
  [[ -s "$fasta" ]] || { echo "[deeparg] skip $context/$label (empty/missing)"; return 0; }

  local ctx_dir="$OUT/$context"; mkdir -p "$ctx_dir"
  local base; base=$(basename "$fasta" | sed 's/\.\(fasta\|fa\|fna\)$//')
  local outpfx="$ctx_dir/${label}"

  # skip if already done for this context
  [[ -s "${outpfx}.mapping.ARG" && "$FORCE" != "1" ]] && {
    echo "[deeparg] skip $context/$label (cached)"; return 0; }

  echo "[deeparg] scanning $context/$label: $fasta"
  deeparg predict \
    --model LS \
    -i "$fasta" \
    -o "$outpfx" \
    -d "$DEEPARG_DATA" \
    --arg-alignment-identity 50 \
    --arg-alignment-evalue 1e-10 \
    --arg-num-alignments-per-entry 1000 \
    2>&1 || { echo "[deeparg] WARNING: $context/$label failed"; return 0; }

  if [[ -s "${outpfx}.mapping.ARG" ]]; then
    local n; n=$(tail -n +2 "${outpfx}.mapping.ARG" | wc -l)
    echo "[deeparg] $context/$label: $n ARG predictions"
  fi
}

total=0

# ---- 1) Contigs ----
CONTIGS="$ROOT/SRA/assemblies/spades/$SAMPLE/contigs.fasta"
[[ -s "$CONTIGS" ]] || CONTIGS="$ROOT/SRA/assemblies/metaviralspades/$SAMPLE/contigs.fasta"
if [[ -s "$CONTIGS" ]]; then
  _run_deeparg "$CONTIGS" "contigs" "contigs"
  ((total++)) || true
fi

# ---- 2) MAGs ----
shopt -s nullglob
MAG_BINS=( "$ROOT/SRA/binning/metabat2/$SAMPLE"/bin.*.fa )
# also try magscot if it exists
MAG_BINS+=( "$ROOT/SRA/binning/magscot/$SAMPLE"/bin.*.fa )
shopt -u nullglob
if [[ ${#MAG_BINS[@]} -gt 0 ]]; then
  for bin in "${MAG_BINS[@]}"; do
    [[ -s "$bin" ]] || continue
    label=$(basename "$bin" .fa)
    _run_deeparg "$bin" "mags" "$label"
    ((total++)) || true
  done
else
  echo "[deeparg] skip MAGs (no bins found for $SAMPLE)"
fi

# ---- 3) Plasmids ----
PLASMIDS=""
for p in "$ROOT/SRA/plasmids/$SAMPLE/union/plasmids_union.fasta" \
         "$ROOT/SRA/plasmids/$SAMPLE/union/plasmids_concat.fasta" \
         "$ROOT/SRA/plasmids/$SAMPLE/union/"*.fasta; do
  [[ -s "$p" ]] && { PLASMIDS="$p"; break; }
done
if [[ -n "$PLASMIDS" ]]; then
  _run_deeparg "$PLASMIDS" "plasmids" "plasmids"
  ((total++)) || true
else
  echo "[deeparg] skip plasmids (no union FASTA for $SAMPLE)"
fi

# ---- 4) Viruses ----
VIRUSES=""
for v in "$ROOT/SRA/viruses/$SAMPLE/union/viruses_union.fasta" \
         "$ROOT/SRA/viruses/$SAMPLE/union/viruses_concat.fasta" \
         "$ROOT/SRA/viruses/$SAMPLE/union/"*.fasta; do
  [[ -s "$v" ]] && { VIRUSES="$v"; break; }
done
if [[ -n "$VIRUSES" ]]; then
  _run_deeparg "$VIRUSES" "viruses" "viruses"
  ((total++)) || true
else
  echo "[deeparg] skip viruses (no union FASTA for $SAMPLE)"
fi

echo "[deeparg] scanned $total context(s)"
date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
echo "[deeparg] done -> $OUT"
