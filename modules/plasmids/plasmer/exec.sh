#!/usr/bin/env bash
# modules/plasmids/plasmer/exec.sh
# Plasmer (k-mer + genomic-feature plasmid predictor) on the metaSPAdes contigs.
# Output: SRA/plasmids/<s>/plasmer/results/<s>.plasmer.predPlasmids.fa (+ predClass/
# predProb tables); picked up by union_cluster with prefix "plasmer".
# Plasmer loads its k-mer DB into RAM (>= 32 GB): on smaller machines the step is
# skipped with a message (PLASMER_FORCE=1 / PLASMER_MIN_MEM_GB override).
# DB: refdata/plasmer (plasmerMainDB + customizedKraken2DB, Zenodo 7030675,
# ~14 GB packed), fetched on first use.
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"
OUT="$ROOT/SRA/plasmids/$SAMPLE/plasmer"; mkdir -p "$OUT"
DONE="$OUT/.done"
[[ "$FORCE" == "1" ]] && rm -f "$DONE"
[[ -s "$DONE" ]] && { echo "[plasmer] cached"; exit 0; }

IN="$ROOT/SRA/assemblies/spades/$SAMPLE/contigs.fasta"
[[ -s "$IN" ]] || { echo "[plasmer] no input: $IN"; exit 1; }
command -v Plasmer >/dev/null 2>&1 || { echo "[plasmer] Plasmer not in PATH"; exit 3; }

# shellcheck source=bin/lib/cfg.sh
source "$ROOT/bin/lib/cfg.sh"
MINLEN="$(cfg modules.plasmids.plasmer_min_len 500)"

DB="${PLASMER_DB:-$ROOT/refdata/plasmer}"
if [[ ! -f "$DB/.extracted" ]]; then
  echo "[plasmer] DB missing -> downloading into $DB (~14 GB packed)"
  mkdir -p "$DB/tmp"
  for f in plasmerMainDB.tar.xz customizedKraken2DB.tar.xz; do
    curl -fL -C - -o "$DB/tmp/$f" "https://zenodo.org/records/7030675/files/$f?download=1"
    tar -xJf "$DB/tmp/$f" -C "$DB"
    rm -f "$DB/tmp/$f"
  done
  rmdir "$DB/tmp" 2>/dev/null || true
  date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DB/.extracted"
fi

echo "[plasmer] Plasmer -g $IN -d $DB -t $CPUS -m $MINLEN -o $OUT"
Plasmer -g "$IN" -p "$SAMPLE" -d "$DB" -t "$CPUS" -m "$MINLEN" -o "$OUT"

PRED="$OUT/results/${SAMPLE}.plasmer.predPlasmids.fa"
[[ -f "$PRED" ]] || { echo "[plasmer] ERROR: no $PRED"; exit 4; }
echo "[plasmer] $(grep -c '^>' "$PRED" || true) predicted plasmid contigs -> $PRED"
date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
