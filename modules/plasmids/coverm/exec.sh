#!/usr/bin/env bash
# Map trimmed reads to global dereplicated plasmids and quantify with CoverM
set -euo pipefail
ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

FA="$ROOT/SRA/plasmids/_global/cluster/plasmids_derep.global.fasta"
R1=( "$ROOT"/SRA/reads/trimmed/"$SAMPLE"/*_R1*.fastq.gz )
R2=( "$ROOT"/SRA/reads/trimmed/"$SAMPLE"/*_R2*.fastq.gz )
OUT="$ROOT/SRA/plasmids/$SAMPLE/coverm"; mkdir -p "$OUT"

[[ -s "$FA" ]] || { echo "[coverm] global derep fasta missing: $FA"; exit 2; }
[[ -n "${R1[0]:-}" ]] || { echo "[coverm] no reads found"; exit 3; }
[[ -n "${R2[0]:-}" ]] || SINGLE_END=1

# Tunables
ID="${MIN_ID:-0.95}"            # min read percent identity
ALN="${MIN_ALN_FRAC:-0.75}"      # min fraction of each read aligned
MAPPER="${MAPPER:-minimap2}"     # minimap2|bowtie2|bwa-mem
DONE="$OUT/.done.id${ID}.aln${ALN}.${MAPPER}"

# Mapper availability check
ok_mapper=0
case "$MAPPER" in
  minimap2) command -v minimap2 >/dev/null && ok_mapper=1 ;;
  bowtie2)  command -v bowtie2  >/dev/null && ok_mapper=1 ;;
  bwa-mem)  command -v bwa      >/dev/null && ok_mapper=1 ;;
esac
[[ $ok_mapper -eq 1 ]] || { echo "[coverm] mapper '$MAPPER' not available. Install minimap2 or set MAPPER=bowtie2/bwa-mem."; exit 4; }

# Rerun logic
if [[ "$FORCE" != "1" && -s "$DONE" && -s "$OUT/coverm.tsv" ]]; then
  echo "[coverm] cached -> $OUT/coverm.tsv"; exit 0
fi

# Build input arguments
args=( contig --reference "$FA"
       --methods covered_fraction tpm rpkm
       --threads "$CPUS"
       --mapper "$MAPPER"
       --min-read-percent-identity "$ID"
       --min-read-aligned-percent "$ALN"
       --exclude-supplementary
     )

if [[ "${SINGLE_END:-0}" -eq 1 ]]; then
  # single-end fallback: pass all R1 as separate files
  reads=()
  for f in "${R1[@]}"; do [[ -s "$f" ]] && reads+=( --single "$f" ); done
  [[ ${#reads[@]} -gt 0 ]] || { echo "[coverm] no valid SE reads"; exit 3; }
  args+=( "${reads[@]}" )
else
  # paired-end: pair by index; support multi-lane
  [[ ${#R1[@]} -eq ${#R2[@]} ]] || { echo "[coverm] R1/R2 count mismatch"; exit 3; }
  for i in "${!R1[@]}"; do
    [[ -s "${R1[$i]}" && -s "${R2[$i]}" ]] && args+=( --coupled "${R1[$i]}" "${R2[$i]}" )
  done
fi

# Output
args+=( --output-file "$OUT/coverm.tsv" )

echo "[coverm] running: coverm ${args[*]}"
coverm "${args[@]}"

date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
echo "[coverm] done -> $OUT/coverm.tsv"
