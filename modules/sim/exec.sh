#!/usr/bin/env bash
set -euo pipefail
: "${ROOT:?ROOT not set}"
: "${SET_NAME:?SET_NAME not set}"
: "${N_READS:?N_READS not set}"
: "${READ_LEN:?READ_LEN not set}"
: "${PLAN:?PLAN not set}"     # tables from scripts/sim_plan.py
: "${MODEL:=HiSeq}"           # kept for compatibility, unused by ART
: "${CPUS:=8}"                # unused by ART
: "${SEED:=42}"               # ART base seed (per sample and reference: SEED + offset)

# env helpers
# shellcheck source=bin/lib/env.sh
source "$ROOT/bin/lib/env.sh"
ensure_env_by_module "sim"   # sim_env must have: art, ncbi-datasets-cli, entrez-direct, unzip, yq, pyyaml

REFDIR="$ROOT/refdata/testsets/$SET_NAME"
mkdir -p "$REFDIR/refs/hosts" "$REFDIR/refs/plasmids" "$REFDIR/refs/phages" "$REFDIR/sim" "$REFDIR/truth"

dl() {
  acc="$1"; out="$2"
  [[ -s "$out/$acc.fna" ]] && return
  mkdir -p "$out"

  if [[ "$acc" =~ ^GCF_|^GCA_ ]]; then
    run_in_env sim_env datasets download genome accession "$acc" --include genome --filename "$out/$acc.zip"
    unzip -qo "$out/$acc.zip" -d "$out/_unz"
    fna=$(find "$out/_unz" -type f \( -name "*genomic.fna" -o -name "*.fna" -o -name "*.fa" \) | head -n1) || true
    [[ -n "${fna:-}" ]] || { echo "[error] no FASTA inside genome package for $acc"; exit 1; }
    cp "$fna" "$out/$acc.fna"
    rm -rf "$out/_unz" "$out/$acc.zip"
  else
    run_in_env sim_env bash -lc "efetch -db nuccore -id '$acc' -format fasta > '$out/$acc.fna'"
    [[ -s "$out/$acc.fna" ]] || { echo "[error] efetch failed for $acc"; exit 1; }
  fi
}

# downloads (one folder per reference type)
while IFS=$'\t' read -r acc type _host; do
  dl "$acc" "$REFDIR/refs/${type}s"
done < <(tail -n +2 "$PLAN/references.tsv")

# combine fasta and write the truth table per sequence record. Records keep the first word
# of their FASTA header as name, so alignments map back to them. A record inside a genome
# assembly (GCF_/GCA_) whose header says "plasmid" is typed plasmid with the assembly as host.
COMBINED="$REFDIR/refs/combined.fna"; : > "$COMBINED"
TRUTH="$REFDIR/truth/provenance.tsv"
printf "type\taccession\tset\thost\tsource\n" > "$TRUTH"
while IFS=$'\t' read -r acc type host; do
  awk -v acc="$acc" -v type="$type" -v host="$host" -v set="$SET_NAME" -v truth="$TRUTH" '
    function flush(){ if (name != "") { print ">" name; for (i = 1; i <= length(seq); i += 70) print substr(seq, i, 70) } }
    /^>/ { flush(); hdr = substr($0, 2); split(hdr, w, " "); name = w[1]; seq = "";
           t = type; h = host
           if (acc ~ /^GC[AF]_/ && tolower(hdr) ~ /plasmid/) { t = "plasmid"; h = acc }
           printf "%s\t%s\t%s\t%s\t%s\n", t, name, set, h, acc >> truth
           next }
    { seq = seq $0 }
    END { flush() }' "$REFDIR/refs/${type}s/$acc.fna" >> "$COMBINED"
done < <(tail -n +2 "$PLAN/references.tsv")

# per-sample weights, per configured accession (all records of an assembly share its weight)
cp "$PLAN/weights.tsv" "$REFDIR/truth/weights.tsv"

# ensure ART available
run_in_env sim_env art_illumina --help >/dev/null 2>&1 || run_in_env sim_env micromamba install -y -n sim_env -c bioconda art

# one read set per sample. A single-sample config keeps the historical layout
# SRA/reads/<set>/metagenome_R{1,2}.fastq.gz; every sample of a run belongs to group <set>.
MANIFEST="$REFDIR/sim/manifest.tsv"
printf "sample\tfq1\tfq2\tgroup\n" > "$MANIFEST"
sidx=0
while read -r sample; do
  sidx=$((sidx + 1))
  OUT="$REFDIR/sim/$sample/metagenome"; mkdir -p "$(dirname "$OUT")"
  ABUND_FILE="$REFDIR/sim/$sample/abundances.tsv"
  awk -F'\t' -v s="$sample" 'NR>1 && $1==s {print $2"\t"$3}' "$PLAN/weights.tsv" > "$ABUND_FILE"
  export REFDIR ABUND_FILE N_READS READ_LEN SEED OUT
  export SEED_OFFSET=$(( (sidx - 1) * 1000 ))

  # simulate reads with ART per reference based on weights
  run_in_env sim_env bash -lc '
set -euo pipefail
rm -f "${OUT}_R1.fastq" "${OUT}_R2.fastq" || true
idx=0
while read -r acc w; do
  idx=$((idx + 1))
  [[ -z "${w:-}" ]] && continue
  awk -v w="$w" "BEGIN{exit !(w+0 > 0)}" || continue

  # locate fasta
  f="$(ls "$REFDIR"/refs/*/"${acc}.fna" 2>/dev/null | head -n1 || true)"
  [[ -s "$f" ]] || { echo "[err] fasta for $acc not found"; exit 1; }

  # genome size
  gs=$(awk "/^>/{next}{n+=length} END{print n}" "$f")
  [[ "$gs" -gt 0 ]] || { echo "[err] zero length for $acc"; exit 1; }

  # per-ref fold to match N_READS * 2 * READ_LEN * weight / gs
  fold=$(awk -v n="$N_READS" -v L="$READ_LEN" -v w="$w" -v gs="$gs" \
         "BEGIN{printf \"%.6f\", (n*2*L*w)/gs}")

  # fixed seed per sample and reference so the same config always yields identical reads
  art_illumina -ss HS25 -na -p -l "$READ_LEN" -f "$fold" -m 300 -s 30 -rs $((SEED + SEED_OFFSET + idx)) \
    -i "$f" -o "$OUT.$acc." >/dev/null
done < "$ABUND_FILE"

# concatenate and gzip
shopt -s nullglob
r1_parts=( "$OUT".*.1.fq )
r2_parts=( "$OUT".*.2.fq )
[[ ${#r1_parts[@]} -gt 0 ]] || { echo "[err] no R1 part files"; exit 1; }
[[ ${#r2_parts[@]} -gt 0 ]] || { echo "[err] no R2 part files"; exit 1; }

cat "${r1_parts[@]}" > "${OUT}_R1.fastq"
cat "${r2_parts[@]}" > "${OUT}_R2.fastq"

if command -v pigz >/dev/null 2>&1; then
  pigz -f "${OUT}_R1.fastq" "${OUT}_R2.fastq"
else
  gzip -f "${OUT}_R1.fastq" "${OUT}_R2.fastq"
fi

# optional cleanup of per-ref chunks
rm -f "$OUT".*.1.fq "$OUT".*.2.fq || true
'

  # link outputs
  OUTREADS="$ROOT/SRA/reads/$sample"; mkdir -p "$OUTREADS"
  ln -sf "${OUT}_R1.fastq.gz" "$OUTREADS/metagenome_R1.fastq.gz"
  ln -sf "${OUT}_R2.fastq.gz" "$OUTREADS/metagenome_R2.fastq.gz"
  printf "%s\t%s\t%s\t%s\n" "$sample" "${OUT}_R1.fastq.gz" "${OUT}_R2.fastq.gz" "$SET_NAME" >> "$MANIFEST"
  echo "[sim] $sample -> $OUTREADS/metagenome_R{1,2}.fastq.gz"
done < <(tail -n +2 "$PLAN/weights.tsv" | cut -f1 | awk '!seen[$0]++')

echo "[sim] done: $sidx sample(s), manifest $MANIFEST"
