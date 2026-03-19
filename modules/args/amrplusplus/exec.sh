#!/usr/bin/env bash
# modules/args/amrplusplus/exec.sh
# AMR++ style read-based alignment to MEGARes v3 database
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

# --- collect trimmed reads from qc/fastp ---
shopt -s nullglob
R1=( "$ROOT/SRA/qc/fastp/$SAMPLE/"*trimmed_1*.fastq.gz )
R2=( "$ROOT/SRA/qc/fastp/$SAMPLE/"*trimmed_2*.fastq.gz )
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
COUNTS="$OUT/resistome_counts.tsv"
echo "[amrplusplus] generating resistome count matrix"
samtools idxstats "$BAM" | \
  awk -F'\t' 'BEGIN{OFS="\t"; print "gene","length","mapped_reads","unmapped_reads"} $3>0{print $1,$2,$3,$4}' \
  > "$COUNTS"

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
