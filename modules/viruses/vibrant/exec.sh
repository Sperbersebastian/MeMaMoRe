#!/usr/bin/env bash
# modules/viruses/vibrant/exec.sh
# Run VIBRANT on metaSPAdes contigs
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

# prefer metaSPAdes contigs
IN="$ROOT/SRA/assemblies/spades/$SAMPLE/contigs.fasta"
[[ -s "$IN" ]] || IN="$ROOT/SRA/assemblies/metaviralspades/$SAMPLE/contigs.fasta"
[[ -s "$IN" ]] || { echo "[vibrant] no contigs found for $SAMPLE"; exit 1; }

OUT="$ROOT/SRA/viruses/$SAMPLE/vibrant"; mkdir -p "$OUT"
DONE="$OUT/.done"
[[ "$FORCE" == "1" ]] && rm -f "$DONE"
[[ -s "$DONE" ]] && { echo "[vibrant] cached"; exit 0; }

command -v VIBRANT_run.py >/dev/null 2>&1 || { echo "[vibrant] VIBRANT_run.py not in PATH"; exit 3; }

# VIBRANT requires its HMM database to be setup before first run.
# The databases live on the Box drive (symlinked from the conda env's share dir).
VIBRANT_DB="/media/Box/MeMaMoRe/refdata/vibrant/db"
export VIBRANT_DATA_PATH="$VIBRANT_DB"
if [[ ! -s "$VIBRANT_DB/databases/Pfam-A_v32.HMM.h3i" ]]; then
  echo "[vibrant] DB missing or incomplete -> downloading and setting up"
  command -v download-db.sh >/dev/null 2>&1 || { echo "[vibrant] download-db.sh not found"; exit 4; }
  download-db.sh "$VIBRANT_DB" || { echo "[vibrant] DB setup failed"; exit 5; }
fi

# VIBRANT expects the input basename for output naming
IN_BASE="$(basename "$IN" .fasta)"

cmd=( VIBRANT_run.py
  -i "$IN"
  -folder "$OUT"
  -t "$CPUS"
  -l 1000
)
echo "[vibrant] ${cmd[*]}"
"${cmd[@]}"

# normalise outputs – VIBRANT creates a nested structure
VIBRANT_DIR="$OUT/VIBRANT_${IN_BASE}"
PHAGE_FA="$VIBRANT_DIR/VIBRANT_phages_${IN_BASE}/${IN_BASE}.phages_combined.fna"
if [[ -s "$PHAGE_FA" ]]; then
  ln -sf "$PHAGE_FA" "$OUT/viruses.fasta"
else
  # try alternative output path
  ALT_FA="$(find "$OUT" -name '*.phages_combined.fna' -print -quit 2>/dev/null || true)"
  if [[ -s "${ALT_FA:-}" ]]; then
    ln -sf "$ALT_FA" "$OUT/viruses.fasta"
  else
    echo "[vibrant] WARNING: no viral sequences found"
  fi
fi

date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
echo "[vibrant] done -> $OUT"
