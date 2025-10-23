#!/usr/bin/env bash
# modules/plasmids/genomad/exec.sh
# Run geNomad end-to-end on contigs, auto-download DB if missing.
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

# ---------- inputs ----------
IN_SPADES="$ROOT/SRA/assemblies/spades/$SAMPLE/contigs.fasta"
IN_MPS="$ROOT/SRA/assemblies/metaplasmidspades/$SAMPLE/contigs.fasta"
if   [[ -s "$IN_SPADES" ]]; then IN="$IN_SPADES"
elif [[ -s "$IN_MPS"   ]]; then IN="$IN_MPS"
else
  echo "[genomad] no contigs: $IN_SPADES or $IN_MPS"; exit 1
fi

# DB path can be overridden via GENOMAD_DB
DB="${GENOMAD_DB:-$ROOT/refdata/genomad/genomad_db}"
DB_BASE="$(dirname "$DB")"

OUT="$ROOT/SRA/plasmids/$SAMPLE/genomad"
mkdir -p "$OUT" "$DB_BASE"

command -v genomad >/dev/null 2>&1 || { echo "[genomad] genomad not in PATH"; exit 3; }

# ---------- ensure database ----------
if [[ ! -f "$DB/version.txt" ]]; then
  echo "[genomad] DB missing -> downloading to $DB_BASE"
  ( cd "$DB_BASE" && genomad download-database . )
  [[ -f "$DB/version.txt" ]] || { echo "[genomad] DB download failed at $DB"; exit 2; }
  echo "[genomad] DB ready: $DB"
fi

# ---------- caching ----------
DONE="$OUT/.done"
[[ "$FORCE" == "1" ]] && rm -f "$DONE"
if [[ -s "$DONE" ]]; then
  echo "[genomad] cached"; exit 0
fi

# ---------- run ----------
cmd=( genomad end-to-end "$IN" "$OUT" "$DB" --threads "$CPUS" )
echo "[genomad] ${cmd[*]}"
"${cmd[@]}"

# ---------- normalize outputs ----------
# TSVs
[[ -s "$OUT/genomad_plasmid_prediction.tsv" && ! -e "$OUT/plasmids.tsv" ]] && \
  ln -sf "genomad_plasmid_prediction.tsv" "$OUT/plasmids.tsv"
[[ -s "$OUT/genomad_virus_prediction.tsv" && ! -e "$OUT/viruses.tsv" ]] && \
  ln -sf "genomad_virus_prediction.tsv" "$OUT/viruses.tsv"

# FASTA
[[ -s "$OUT/genomad_plasmids.fna" && ! -e "$OUT/plasmids.fna" ]] && \
  ln -sf "genomad_plasmids.fna" "$OUT/plasmids.fna"
[[ -s "$OUT/genomad_viruses.fna" && ! -e "$OUT/viruses.fna" ]] && \
  ln -sf "genomad_viruses.fna" "$OUT/viruses.fna"

date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
echo "[genomad] done -> $OUT"
