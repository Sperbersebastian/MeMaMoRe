#!/usr/bin/env bash
# modules/viruses/checkv/exec.sh
# Run CheckV on deduplicated virus union set
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

IN="$ROOT/SRA/viruses/$SAMPLE/union/viruses_union.fasta"
OUT="$ROOT/SRA/viruses/$SAMPLE/checkv"; mkdir -p "$OUT"

[[ -s "$IN" ]] || { echo "[checkv] no input: $IN"; exit 1; }

DONE="$OUT/.done"
[[ "$FORCE" == "1" ]] && rm -f "$DONE"
[[ -s "$DONE" ]] && { echo "[checkv] cached"; exit 0; }

command -v checkv >/dev/null 2>&1 || { echo "[checkv] checkv not in PATH"; exit 3; }

# CheckV database (auto-download if missing)
CHECKV_DB="${CHECKVDB:-${CHECKV_DB:-$ROOT/refdata/checkv/checkv-db-v1.5}}"
export CHECKVDB="$CHECKV_DB"

if [[ ! -s "$CHECKV_DB/genome_db/checkv_reps.faa" ]]; then
  echo "[checkv] DB missing -> downloading via curl"
  mkdir -p "$(dirname "$CHECKV_DB")"
  curl -C - -L -o checkv-db-v1.5.tar.gz https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz
  tar -xzf checkv-db-v1.5.tar.gz -C "$(dirname "$CHECKV_DB")"
  rm -f checkv-db-v1.5.tar.gz
fi

# CheckV DB requires running diamond makedb on the .faa file
if [[ -s "$CHECKV_DB/genome_db/checkv_reps.faa" && ! -s "$CHECKV_DB/genome_db/checkv_reps.dmnd" ]]; then
  echo "[checkv] Building diamond database"
  diamond makedb --in "$CHECKV_DB/genome_db/checkv_reps.faa" -d "$CHECKV_DB/genome_db/checkv_reps"
fi

[[ -d "$CHECKV_DB" && -s "$CHECKV_DB/genome_db/checkv_reps.dmnd" ]] || { echo "[checkv] DB extraction or diamond build failed"; exit 2; }
export CHECKVDB="$CHECKV_DB"

cmd=( checkv end_to_end "$IN" "$OUT" -t "$CPUS" -d "$CHECKV_DB" )
echo "[checkv] ${cmd[*]}"
"${cmd[@]}"

# Summary
if [[ -s "$OUT/quality_summary.tsv" ]]; then
  echo "[checkv] quality summary:"
  head -n 5 "$OUT/quality_summary.tsv"
fi

date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
echo "[checkv] done -> $OUT"
