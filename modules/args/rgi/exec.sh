#!/usr/bin/env bash
# modules/args/rgi/exec.sh
# Run RGI (CARD) on multiple genomic contexts: contigs, MAGs, plasmids, viruses
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

OUT="$ROOT/SRA/args/$SAMPLE/rgi"; mkdir -p "$OUT"
DONE="$OUT/.done"
[[ "$FORCE" == "1" ]] && rm -f "$DONE"
[[ -s "$DONE" ]] && { echo "[rgi] cached"; exit 0; }

command -v rgi >/dev/null 2>&1 || { echo "[rgi] rgi not in PATH"; exit 3; }

# Auto-setup CARD database on first run
# RGI --local stores/looks for "localDB" in the current working directory
CARD_DATA="${CARD_DATA:-$ROOT/refdata/card}"
mkdir -p "$CARD_DATA"
cd "$CARD_DATA"

if [[ ! -d "$CARD_DATA/localDB" ]]; then
  echo "[rgi] CARD local database not found -> downloading and loading"
  CARD_JSON="$CARD_DATA/card.json"
  if [[ ! -s "$CARD_JSON" ]]; then
    echo "[rgi] downloading CARD data..."
    curl -sSL "https://card.mcmaster.ca/latest/data" -o "$CARD_DATA/card-data.tar.bz2" \
      || { echo "[rgi] CARD download failed"; exit 5; }
    tar xjf "$CARD_DATA/card-data.tar.bz2" -C "$CARD_DATA" \
      || { echo "[rgi] CARD extraction failed"; exit 5; }
    # card.json may be nested — find it
    if [[ ! -s "$CARD_JSON" ]]; then
      found=$(find "$CARD_DATA" -name "card.json" -type f | head -1)
      [[ -n "$found" ]] && ln -sf "$found" "$CARD_JSON"
    fi
  fi
  [[ -s "$CARD_JSON" ]] || { echo "[rgi] card.json not found after extraction"; exit 5; }
  rgi load --card_json "$CARD_JSON" --local || { echo "[rgi] CARD load failed"; exit 5; }
  echo "[rgi] CARD database loaded"
fi

NCPUS=$(( CPUS < $(nproc) ? CPUS : $(nproc) ))

# ---- helper: run RGI on one FASTA ----
_run_rgi() {
  local fasta="$1" context="$2" label="$3"
  [[ -s "$fasta" ]] || { echo "[rgi] skip $context/$label (empty/missing)"; return 0; }

  local ctx_dir="$OUT/$context"; mkdir -p "$ctx_dir"
  local outpfx="$ctx_dir/${label}"

  # skip if already done for this context
  [[ -s "${outpfx}.txt" && "$FORCE" != "1" ]] && {
    echo "[rgi] skip $context/$label (cached)"; return 0; }

  echo "[rgi] scanning $context/$label: $fasta"
  # run from CARD_DATA dir so --local finds localDB
  rgi main \
    -i "$fasta" \
    -o "$outpfx" \
    -t contig \
    -a DIAMOND \
    -n "$NCPUS" \
    --include_loose \
    --local \
    --clean \
    2>&1 || { echo "[rgi] WARNING: $context/$label failed"; return 0; }

  if [[ -s "${outpfx}.txt" ]]; then
    local n; n=$(tail -n +2 "${outpfx}.txt" | wc -l)
    echo "[rgi] $context/$label: $n AMR gene hits"
  fi
}

total=0

# ---- 1) Contigs ----
CONTIGS="$ROOT/SRA/assemblies/spades/$SAMPLE/contigs.fasta"
[[ -s "$CONTIGS" ]] || CONTIGS="$ROOT/SRA/assemblies/metaviralspades/$SAMPLE/contigs.fasta"
if [[ -s "$CONTIGS" ]]; then
  _run_rgi "$CONTIGS" "contigs" "contigs"
  ((total++)) || true
fi

# ---- 2) MAGs ----
shopt -s nullglob
# refined MAGScoT bins if present, else raw MetaBAT2 bins
MAG_BINS=( "$ROOT/SRA/binning/magscot/$SAMPLE/bins"/*.fa )
[[ ${#MAG_BINS[@]} -gt 0 ]] || MAG_BINS=( "$ROOT/SRA/binning/metabat2/$SAMPLE"/bin.*.fa )
shopt -u nullglob
if [[ ${#MAG_BINS[@]} -gt 0 ]]; then
  for bin in "${MAG_BINS[@]}"; do
    [[ -s "$bin" ]] || continue
    label=$(basename "$bin" .fa)
    _run_rgi "$bin" "mags" "$label"
    ((total++)) || true
  done
else
  echo "[rgi] skip MAGs (no bins found for $SAMPLE)"
fi

# ---- 3) Plasmids ----
PLASMIDS=""
for p in "$ROOT/SRA/plasmids/$SAMPLE/cluster/plasmids_derep.$SAMPLE.fasta" \
         "$ROOT/SRA/plasmids/$SAMPLE/union/plasmids_union.fasta" \
         "$ROOT/SRA/plasmids/$SAMPLE/union/plasmids_concat.fasta" \
         "$ROOT/SRA/plasmids/$SAMPLE/union/"*.fasta; do
  [[ -s "$p" ]] && { PLASMIDS="$p"; break; }
done
if [[ -n "$PLASMIDS" ]]; then
  _run_rgi "$PLASMIDS" "plasmids" "plasmids"
  ((total++)) || true
else
  echo "[rgi] skip plasmids (no union FASTA for $SAMPLE)"
fi

# ---- 4) Viruses ----
VIRUSES=""
for v in "$ROOT/SRA/viruses/$SAMPLE/union/viruses_union.fasta" \
         "$ROOT/SRA/viruses/$SAMPLE/union/viruses_concat.fasta" \
         "$ROOT/SRA/viruses/$SAMPLE/union/"*.fasta; do
  [[ -s "$v" ]] && { VIRUSES="$v"; break; }
done
if [[ -n "$VIRUSES" ]]; then
  _run_rgi "$VIRUSES" "viruses" "viruses"
  ((total++)) || true
else
  echo "[rgi] skip viruses (no union FASTA for $SAMPLE)"
fi

echo "[rgi] scanned $total context(s)"
date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
echo "[rgi] done -> $OUT"
