#!/usr/bin/env bash
# modules/args/rgi/exec.sh
# Run RGI (CARD, protein mode) once on the central gene calls (genes.faa)
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

# ---- input: central gene calls; RGI protein mode on genes.faa ----
# (protein mode cannot evaluate rRNA mutation models; all others apply)
GENES="$ROOT/SRA/assemblies/genes/$SAMPLE/genes.faa"
[[ -s "$GENES" ]] || { echo "[rgi] missing $GENES (run: bin/main.sh run assembly --sample $SAMPLE)"; exit 2; }

# Every gene is scanned once; context is assigned later from the contig ID.
CTX_DIR="$OUT/genes"; mkdir -p "$CTX_DIR"
OUTPFX="$CTX_DIR/genes"
if [[ -s "${OUTPFX}.txt" && "$FORCE" != "1" && "${OUTPFX}.txt" -nt "$GENES" ]]; then
  echo "[rgi] cached -> ${OUTPFX}.txt"
else
  # RGI rejects '*' stop symbols in protein input
  sed '/^>/!s/\*//g' "$GENES" > "$CTX_DIR/genes.nostop.faa"
  echo "[rgi] scanning $(grep -c '^>' "$GENES") proteins: $GENES"
  rgi main \
    -i "$CTX_DIR/genes.nostop.faa" \
    -o "$OUTPFX" \
    -t protein \
    -a DIAMOND \
    -n "$NCPUS" \
    --include_loose \
    --local \
    --clean \
    2>&1
  rm -f "$CTX_DIR/genes.nostop.faa"
  [[ -f "${OUTPFX}.txt" ]] || { echo "[rgi] ERROR: no ${OUTPFX}.txt"; exit 4; }
fi
echo "[rgi] $(tail -n +2 "${OUTPFX}.txt" | wc -l) AMR gene hits (incl. Loose)"

date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
echo "[rgi] done -> $OUT"
