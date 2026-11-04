#!/usr/bin/env bash
# modules/args/amrplusplus/exec.sh
# AMR++ style read-based alignment to MEGARes v3 database
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

# --- collect trimmed reads from qc/fastp ---
shopt -s nullglob
# shellcheck source=bin/lib/reads.sh
source "$ROOT/bin/lib/reads.sh"
# fastp output, or host-cleaned reads when modules.qc.host_removal = human
_rr="$(qc_reads "$SAMPLE")" || exit 2
IFS=$'\t' read -r _r1 _r2 <<<"$_rr"
R1=(); R2=()
[[ -s "$_r1" ]] && R1=( "$_r1" )
[[ -n "$_r2" && -s "$_r2" ]] && R2=( "$_r2" )
shopt -u nullglob
[[ ${#R1[@]} -gt 0 ]] || { echo "[amrplusplus] no QC'd reads for $SAMPLE"; exit 1; }

OUT="$ROOT/SRA/args/$SAMPLE/amrplusplus"; mkdir -p "$OUT"
DONE="$OUT/.done"
[[ "$FORCE" == "1" ]] && rm -f "$DONE"
[[ -s "$DONE" ]] && { echo "[amrplusplus] cached"; exit 0; }

command -v bwa >/dev/null 2>&1 || { echo "[amrplusplus] bwa not in PATH"; exit 3; }
command -v samtools >/dev/null 2>&1 || { echo "[amrplusplus] samtools not in PATH"; exit 3; }

# --- MEGARes v3 database ---
MEGARES_DIR="${MEGARES_DIR:-$ROOT/refdata/amrplusplus/megares_v3}"
MEGARES_FA="$MEGARES_DIR/megares_database_v3.00.fasta"
MEGARES_ANN="$MEGARES_DIR/megares_annotations_v3.00.csv"

if [[ ! -s "$MEGARES_FA" ]]; then
  echo "[amrplusplus] downloading MEGARes v3 database"
  mkdir -p "$MEGARES_DIR"
  curl -sSL "https://www.meglab.org/downloads/megares_v3.00/megares_database_v3.00.fasta" \
    -o "$MEGARES_FA" || { echo "[amrplusplus] DB download failed"; exit 5; }
  curl -sSL "https://www.meglab.org/downloads/megares_v3.00/megares_annotations_v3.00.csv" \
    -o "$MEGARES_ANN" || true
fi

# --- BWA index ---
if [[ ! -s "${MEGARES_FA}.bwt" ]]; then
  echo "[amrplusplus] indexing MEGARes database"
  bwa index "$MEGARES_FA"
fi

# --- Align reads ---
BAM="$OUT/megares_aligned.bam"
echo "[amrplusplus] aligning reads to MEGARes..."
if [[ ${#R2[@]} -gt 0 ]]; then
  bwa mem -t "$CPUS" "$MEGARES_FA" "${R1[0]}" "${R2[0]}" 2>/dev/null | \
    samtools view -bS -F 4 -@ "$CPUS" | \
    samtools sort -@ "$CPUS" -o "$BAM"
else
  bwa mem -t "$CPUS" "$MEGARES_FA" "${R1[0]}" 2>/dev/null | \
    samtools view -bS -F 4 -@ "$CPUS" | \
    samtools sort -@ "$CPUS" -o "$BAM"
fi
samtools index "$BAM"

# --- Generate resistome counts ---
# All genes with >= 1 read go to resistome_counts.all.tsv. Like AMR++ itself,
# resistome_counts.tsv keeps only genes covered over >= AMRPP_MIN_GF % of their
# length (default 80) and drops MEGARes entries that need SNP confirmation.
COUNTS="$OUT/resistome_counts.tsv"
ALLCOUNTS="$OUT/resistome_counts.all.tsv"
MIN_GF="${AMRPP_MIN_GF:-80}"
echo "[amrplusplus] generating resistome count matrix (gene fraction >= ${MIN_GF}%)"
samtools coverage "$BAM" | \
  awk -F'\t' 'BEGIN{OFS="\t"; print "gene","length","mapped_reads","gene_fraction"}
               !/^#/ && $4>0 {print $1, $3-$2+1, $4, $6}' > "$ALLCOUNTS"
awk -F'\t' -v min="$MIN_GF" 'BEGIN{OFS="\t"} NR==1 || ($4>=min && $1 !~ /RequiresSNPConfirmation/)' \
  "$ALLCOUNTS" > "$COUNTS"
echo "[amrplusplus] $(( $(wc -l < "$ALLCOUNTS") - 1 )) genes with reads, $(( $(wc -l < "$COUNTS") - 1 )) pass the gene-fraction filter"

NGENES=$(tail -n +2 "$COUNTS" | wc -l)
NREADS=$(tail -n +2 "$COUNTS" | awk '{s+=$3}END{print s+0}')
echo "[amrplusplus] $NGENES genes hit, $NREADS aligned reads"

# --- Annotate with MEGARes classes if annotation file available ---
if [[ -s "$MEGARES_ANN" ]]; then
  ANNOTATED="$OUT/resistome_annotated.tsv"
  python3 -c "
import csv, sys
ann = {}
with open('$MEGARES_ANN') as f:
    reader = csv.reader(f)
    hdr = next(reader)
    for row in reader:
        if len(row) >= 4:
            ann[row[0]] = {'class': row[1], 'mechanism': row[2], 'group': row[3]}

print('gene\tlength\tmapped_reads\tclass\tmechanism\tgroup')
with open('$COUNTS') as f:
    next(f)  # skip header
    for line in f:
        parts = line.strip().split('\t')
        gene = parts[0]
        a = ann.get(gene, {'class':'unknown','mechanism':'unknown','group':'unknown'})
        print(f'{parts[0]}\t{parts[1]}\t{parts[2]}\t{a[\"class\"]}\t{a[\"mechanism\"]}\t{a[\"group\"]}')
" > "$ANNOTATED" 2>/dev/null && echo "[amrplusplus] annotated output -> $ANNOTATED"
fi

date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
echo "[amrplusplus] done -> $OUT"
