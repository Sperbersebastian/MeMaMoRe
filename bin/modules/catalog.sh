#!/usr/bin/env bash
# Project MAG catalog (bin/main.sh run catalog [--samples s1,s2]):
#   1) collect the refined MAGs of all samples (SRA/binning/magscot/<s>/bins/*.fa)
#   2) dRep v4 dereplicate with skani at modules.catalog.ani (0.95), genome
#      quality from CheckM2 via --genomeInfo, quality filter from config
#   3) CoverM genome: abundance of every catalog genome in every sample
#      (one read mapping per sample against the catalog)
#   4) SRA/summary/mag_catalog.tsv (+ SRA/catalog/ working dir)
# Runs inside drep_env (dRep, skani, CoverM).
set -euo pipefail
: "${ROOT:?}"
# shellcheck source=bin/lib/reads.sh
source "$ROOT/bin/lib/reads.sh"

CPUS="$(cfg global.threads 8)"
ANI="$(cfg modules.catalog.ani 0.95)"
COMP="$(cfg modules.catalog.min_completeness 50)"
CONT="$(cfg modules.catalog.max_contamination 10)"
OUT="$ROOT/SRA/catalog"
LOG="$ROOT/logs/catalog.log"
mkdir -p "$OUT" "$ROOT/logs" "$ROOT/SRA/summary"

if [[ -n "${CATALOG_SAMPLES:-}" ]]; then
  IFS=',' read -r -a SAMPLES <<<"$CATALOG_SAMPLES"
else
  mapfile -t SAMPLES < <(tail -n +2 "$ROOT/SRA/reads/manifest.tsv" | cut -f1 | sed '/^$/d')
fi

# ---- 1) genomes + CheckM2 quality (genome,completeness,contamination) ----
G="$OUT/genomes"; rm -rf "$G"; mkdir -p "$G"
INFO="$OUT/genomeInfo.csv"; echo "genome,completeness,contamination" > "$INFO"
for s in "${SAMPLES[@]}"; do
  shopt -s nullglob
  for b in "$ROOT/SRA/binning/magscot/$s/bins"/*.fa; do
    n="${s}__$(basename "$b")"
    ln -sfn "$b" "$G/$n"
  done
  shopt -u nullglob
  q="$ROOT/SRA/binning/checkm2/$s/quality_report.tsv"
  [[ -s "$q" ]] && awk -F'\t' -v s="$s" 'NR>1 {printf "%s__%s.fa,%s,%s\n", s, $1, $2, $3}' "$q" >> "$INFO"
done
NG=$(find "$G" -name '*.fa' | wc -l)
[[ "$NG" -gt 0 ]] || { echo "[catalog] no refined MAGs in ${SAMPLES[*]} (run binning)"; exit 2; }
echo "[catalog] $NG MAGs from ${#SAMPLES[@]} sample(s); dRep skani ANI >= $ANI, completeness >= $COMP, contamination <= $CONT"

# ---- 2) dRep v4 ----
rm -rf "$OUT/drep"
dRep dereplicate "$OUT/drep" -g "$G"/*.fa -p "$CPUS" --S_algorithm skani -sa "$ANI" \
  --genomeInfo "$INFO" -comp "$COMP" -con "$CONT" >>"$LOG" 2>&1
REPS="$OUT/drep/dereplicated_genomes"
NR=$(find "$REPS" -name '*.fa' 2>/dev/null | wc -l)
[[ "$NR" -gt 0 ]] || { echo "[catalog] dRep kept no genome (quality filter?) - see $LOG"; exit 3; }
echo "[catalog] $NR representative genome(s)"

# ---- 3) CoverM abundance per sample ----
AB="$OUT/abundance"; mkdir -p "$AB"
for s in "${SAMPLES[@]}"; do
  rr="$(qc_reads "$s")" || { echo "[catalog] $s: no QC'd reads, no abundance"; continue; }
  IFS=$'\t' read -r R1 R2 <<<"$rr"
  [[ -s "$R1" ]] || { echo "[catalog] $s: no QC'd reads, no abundance"; continue; }
  reads=(-1 "$R1"); [[ -n "$R2" ]] && reads+=(-2 "$R2")
  [[ -n "$R2" ]] || reads=(--single "$R1")
  coverm genome "${reads[@]}" --genome-fasta-directory "$REPS" -x fa \
    --methods relative_abundance mean covered_fraction -t "$CPUS" \
    -o "$AB/$s.tsv" >>"$LOG" 2>&1
done

# ---- 4) table ----
python3 "$ROOT/scripts/mag_catalog.py" "$ROOT" "$OUT" "${SAMPLES[@]}"
