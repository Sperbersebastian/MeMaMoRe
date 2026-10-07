#!/usr/bin/env bash
# modules/viruses/genomad/exec.sh
# Symlink geNomad virus outputs from plasmids run, or run independently
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

OUT="$ROOT/SRA/viruses/$SAMPLE/genomad"; mkdir -p "$OUT"
DONE="$OUT/.done"
[[ "$FORCE" == "1" ]] && rm -f "$DONE"
[[ -s "$DONE" ]] && { echo "[genomad-virus] cached"; exit 0; }

# Try to reuse plasmids geNomad run
PLASMID_GENOMAD="$ROOT/SRA/plasmids/$SAMPLE/genomad"
if [[ -d "$PLASMID_GENOMAD" ]]; then
  echo "[genomad-virus] reusing plasmids geNomad output"
  for f in "$PLASMID_GENOMAD"/*virus*; do
    [[ -e "$f" ]] && ln -sf "$f" "$OUT/$(basename "$f")"
  done
  # Also link the summary output directory if it exists
  # *_summary holds contigs_virus_summary.tsv (read by the GUI)
  for d in "$PLASMID_GENOMAD"/*_find "$PLASMID_GENOMAD"/*_summary; do
    [[ -d "$d" ]] && ln -sfn "$d" "$OUT/$(basename "$d")"
  done
  date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
  echo "[genomad-virus] done (symlinked) -> $OUT"
  exit 0
fi

# Otherwise run geNomad independently on metaSPAdes contigs
IN_SPADES="$ROOT/SRA/assemblies/spades/$SAMPLE/contigs.fasta"
IN_MVS="$ROOT/SRA/assemblies/metaviralspades/$SAMPLE/contigs.fasta"
if   [[ -s "$IN_SPADES" ]]; then IN="$IN_SPADES"
elif [[ -s "$IN_MVS"    ]]; then IN="$IN_MVS"
else
  echo "[genomad-virus] no contigs found"; exit 1
fi

DB="${GENOMAD_DB:-$ROOT/refdata/genomad/genomad_db}"
DB_BASE="$(dirname "$DB")"
mkdir -p "$DB_BASE"

command -v genomad >/dev/null 2>&1 || { echo "[genomad-virus] genomad not in PATH"; exit 3; }

# ensure database
if [[ ! -f "$DB/version.txt" ]]; then
  echo "[genomad-virus] DB missing -> downloading to $DB_BASE"
  ( cd "$DB_BASE" && genomad download-database . )
  [[ -f "$DB/version.txt" ]] || { echo "[genomad-virus] DB download failed"; exit 2; }
fi

cmd=( genomad end-to-end "$IN" "$OUT" "$DB" --threads "$CPUS" )
# GENOMAD_SPLITS>0 lowers MMseqs2 RAM use (needed on machines with <~20 GB RAM)
[[ "${GENOMAD_SPLITS:-0}" -gt 0 ]] && cmd+=( --splits "$GENOMAD_SPLITS" )
echo "[genomad-virus] ${cmd[*]}"
"${cmd[@]}"

# normalise virus outputs – geNomad names files after the input basename
# e.g. contigs_summary/contigs_virus.fna  or  genomad_viruses.fna
VIRUS_FNA=""
for candidate in \
  "$OUT"/genomad_viruses.fna \
  "$OUT"/*_summary/*_virus.fna \
  "$OUT"/*_virus.fna; do
  [[ -s "$candidate" ]] && { VIRUS_FNA="$candidate"; break; }
done
[[ -n "$VIRUS_FNA" && ! -e "$OUT/viruses.fna" ]] && \
  ln -sf "$VIRUS_FNA" "$OUT/viruses.fna"

VIRUS_TSV=""
for candidate in \
  "$OUT"/genomad_virus_prediction.tsv \
  "$OUT"/*_summary/*_virus_summary.tsv; do
  [[ -s "$candidate" ]] && { VIRUS_TSV="$candidate"; break; }
done
[[ -n "$VIRUS_TSV" && ! -e "$OUT/viruses.tsv" ]] && \
  ln -sf "$VIRUS_TSV" "$OUT/viruses.tsv"

date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
echo "[genomad-virus] done -> $OUT"
