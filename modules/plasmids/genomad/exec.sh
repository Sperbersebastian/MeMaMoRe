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
# GENOMAD_SPLITS>0 lowers MMseqs2 RAM use (needed on machines with <~20 GB RAM)
[[ "${GENOMAD_SPLITS:-0}" -gt 0 ]] && cmd+=( --splits "$GENOMAD_SPLITS" )
echo "[genomad] ${cmd[*]}"
"${cmd[@]}"

# ---------- normalize outputs ----------
# geNomad names files after the input basename, e.g. contigs_summary/contigs_plasmid_summary.tsv
_link_first(){ # $1 link name, rest: candidates
  local link="$OUT/$1"; shift
  [[ -e "$link" ]] && return 0
  local c; for c in "$@"; do [[ -s "$c" ]] && { ln -sf "$c" "$link"; return 0; }; done
  return 0
}
_link_first plasmids.tsv "$OUT"/genomad_plasmid_prediction.tsv "$OUT"/*_summary/*_plasmid_summary.tsv
_link_first viruses.tsv  "$OUT"/genomad_virus_prediction.tsv   "$OUT"/*_summary/*_virus_summary.tsv
_link_first plasmids.fna "$OUT"/genomad_plasmids.fna "$OUT"/*_summary/*_plasmid.fna
_link_first viruses.fna  "$OUT"/genomad_viruses.fna  "$OUT"/*_summary/*_virus.fna

date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
echo "[genomad] done -> $OUT"
