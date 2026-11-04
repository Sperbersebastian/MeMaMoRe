#!/usr/bin/env bash
# modules/args/amrfinderplus/exec.sh
# NCBI AMRFinderPlus on the central gene calls (-p genes.faa -g genes.gff, --plus),
# then argNorm maps the hits of AMRFinderPlus, DeepARG and ABRicate/ResFinder to
# CARD ARO terms (RGI reports ARO itself; AMR++ is read-based and not covered).
# Outputs SRA/args/<s>/amrfinderplus/genes.tsv (+ .argnorm.tsv) and
# <tool>/.../*.argnorm.tsv next to the other tools' results.
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"
OUT="$ROOT/SRA/args/$SAMPLE/amrfinderplus"; mkdir -p "$OUT"
DONE="$OUT/.done"
[[ "$FORCE" == "1" ]] && rm -f "$DONE"
[[ -s "$DONE" ]] && { echo "[amrfinderplus] cached"; exit 0; }
command -v amrfinder >/dev/null 2>&1 || { echo "[amrfinderplus] amrfinder not in PATH"; exit 3; }

G="$ROOT/SRA/assemblies/genes/$SAMPLE"
[[ -s "$G/genes.faa" && -s "$G/genes.gff" ]] || { echo "[amrfinderplus] missing $G/genes.{faa,gff} (run assembly)"; exit 2; }

# first run: fetch the AMRFinderPlus DB into the env (~0.2 GB)
amrfinder --database_version >/dev/null 2>&1 || amrfinder -u

# AMRFinderPlus wants bare protein IDs (no '*') and the ID as GFF 'Name='
sed -e '/^>/!s/\*//g' -e 's/^\(>[^ ]*\).*/\1/' "$G/genes.faa" > "$OUT/genes.faa"
awk -F'\t' 'BEGIN{OFS="\t"} /^#/{print; next}
            { if (match($9, /ID=[^;]+/)) $9 = $9 ";Name=" substr($9, RSTART+3, RLENGTH-3); print }' \
  "$G/genes.gff" > "$OUT/genes.gff"

amrfinder -p "$OUT/genes.faa" -g "$OUT/genes.gff" --plus -t "$CPUS" -o "$OUT/genes.tsv"
rm -f "$OUT/genes.faa" "$OUT/genes.gff"
echo "[amrfinderplus] $(awk -F'\t' 'NR>1 && $9=="AMR"' "$OUT/genes.tsv" | wc -l) AMR / $(( $(wc -l < "$OUT/genes.tsv") - 1 )) total hits (--plus)"

# ---- argNorm: ARO terms for the hits of every supported tool ----
if command -v argnorm >/dev/null 2>&1; then
  _norm(){ # tool input output [--db X]
    local tool="$1" in="$2" out="$3"; shift 3
    [[ -s "$in" ]] || return 0
    argnorm "$tool" "$@" -i "$in" -o "$out" 2>&1 | grep -v 'pronto\|_ARO =' || true
  }
  A="$ROOT/SRA/args/$SAMPLE"
  _norm amrfinderplus "$OUT/genes.tsv" "$OUT/genes.argnorm.tsv"
  _norm deeparg "$A/deeparg/genes/genes.mapping.ARG" "$A/deeparg/genes/genes.argnorm.tsv"
  _norm abricate "$A/abricate/resfinder/contigs/contigs.tsv" "$A/abricate/resfinder/contigs/contigs.argnorm.tsv" --db resfinder
else
  echo "[amrfinderplus] argnorm not in PATH -> no ARO normalisation"
fi

date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
echo "[amrfinderplus] done -> $OUT"
