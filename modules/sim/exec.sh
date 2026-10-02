#!/usr/bin/env bash
set -euo pipefail
: "${ROOT:?ROOT not set}"
: "${SET_NAME:?SET_NAME not set}"
: "${N_READS:?N_READS not set}"
: "${READ_LEN:?READ_LEN not set}"
: "${MODEL:?MODEL not set}"   # kept for compatibility, unused by ART
: "${CPUS:=8}"                # unused by ART
: "${SEED:=42}"               # ART base seed (per reference: SEED + index)

# env helpers
# shellcheck source=bin/lib/env.sh
source "$ROOT/bin/lib/env.sh"
ensure_env_by_module "sim"   # sim_env must have: art, ncbi-datasets-cli, entrez-direct, unzip, yq

REFDIR="$ROOT/refdata/testsets/$SET_NAME"
mkdir -p "$REFDIR/refs/hosts" "$REFDIR/refs/plasmids" "$REFDIR/refs/phages" "$REFDIR/sim" "$REFDIR/truth"
OUTREADS="$ROOT/SRA/reads/$SET_NAME"; mkdir -p "$OUTREADS"

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

# downloads
dl "$HOST_HQ" "$REFDIR/refs/hosts"
dl "$HOST_LQ" "$REFDIR/refs/hosts"
for p in ${PLASMIDS}; do dl "$p" "$REFDIR/refs/plasmids"; done
for v in ${PHAGES};   do dl "$v" "$REFDIR/refs/phages";   done

# combine fasta
COMBINED="$REFDIR/refs/combined.fna"; : > "$COMBINED"
for f in "$REFDIR"/refs/hosts/*.fna "$REFDIR"/refs/plasmids/*.fna "$REFDIR"/refs/phages/*.fna; do
  [[ -s "$f" ]] || continue
  awk 'BEGIN{RS=">";FS="\n"} NR>1{
         h=$1; s="";
         for(i=2;i<=NF;i++) s=s $i;
         split(h,a," ");
         print ">" a[1];
         for(i=1;i<=length(s);i+=70) print substr(s,i,70)
       }' "$f" >> "$COMBINED"
done

# abundances
ABUND_FILE="$REFDIR/sim/abundances.tsv"; : > "$ABUND_FILE"
grep "^>" "$COMBINED" | sed 's/^>//' | while read -r acc; do
  var="ABUND_${acc//./_}"
  w="${!var:-0.0}"
  printf "%s\t%s\n" "$acc" "$w" >> "$ABUND_FILE"
done

# export vars for subshell
export REFDIR ABUND_FILE N_READS READ_LEN SEED

# ensure ART available
run_in_env sim_env art_illumina --help >/dev/null 2>&1 || run_in_env sim_env micromamba install -y -n sim_env -c bioconda art

# simulate reads with ART per reference based on weights
run_in_env sim_env bash -lc '
set -euo pipefail
OUT="$REFDIR/sim/metagenome"
rm -f "${OUT}_R1.fastq" "${OUT}_R2.fastq" || true
idx=0

while read -r acc w; do
  [[ -z "${w:-}" || "$w" == "0.0" ]] && continue

  # locate fasta
  f="$REFDIR/refs/${acc}.fna"
  [[ -s "$f" ]] || f="$(ls "$REFDIR"/refs/*/"${acc}.fna" 2>/dev/null | head -n1 || true)"
  [[ -s "$f" ]] || { echo "[err] fasta for $acc not found"; exit 1; }

  # genome size
  gs=$(awk "/^>/{if(n)print n; n=0; next}{n+=length} END{print n}" "$f")
  [[ "$gs" -gt 0 ]] || { echo "[err] zero length for $acc"; exit 1; }

  # per-ref fold to match N_READS * 2 * READ_LEN * weight / gs
  fold=$(awk -v n="$N_READS" -v L="$READ_LEN" -v w="$w" -v gs="$gs" \
         "BEGIN{printf \"%.6f\", (n*2*L*w)/gs}")

  # fixed seed per reference so the same config always yields identical reads
  idx=$((idx + 1))
  art_illumina -ss HS25 -na -p -l '"$READ_LEN"' -f "$fold" -m 300 -s 30 -rs $((SEED + idx)) \
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
ln -sf "$REFDIR/sim/metagenome_R1.fastq.gz" "$OUTREADS/metagenome_R1.fastq.gz"
ln -sf "$REFDIR/sim/metagenome_R2.fastq.gz" "$OUTREADS/metagenome_R2.fastq.gz"

# provenance
{
  printf "type\taccession\tset\n"
  printf "host\t%s\t%s\n" "$HOST_HQ" "$SET_NAME"
  printf "host\t%s\t%s\n" "$HOST_LQ" "$SET_NAME"
  for p in ${PLASMIDS}; do printf "plasmid\t%s\t%s\n" "$p" "$SET_NAME"; done
  for v in ${PHAGES};   do printf "phage\t%s\t%s\n"   "$v" "$SET_NAME"; done
} > "$REFDIR/truth/provenance.tsv"

echo "[sim] done -> $OUTREADS/metagenome_R{1,2}.fastq.gz"
