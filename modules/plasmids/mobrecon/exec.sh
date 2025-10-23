#!/usr/bin/env bash
# modules/plasmids/mobrecon/exec.sh
# Always overwrite safely; works with MOB-suite >=3.1.x
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

IN="$ROOT/SRA/assemblies/spades/$SAMPLE/contigs.fasta"
OUT="$ROOT/SRA/plasmids/$SAMPLE/mobrecon"
DB="${MOBS_DB:-$ROOT/refdata/mobsuite_db}"

[[ -s "$IN" ]] || { echo "[mob_recon] no input: $IN"; exit 1; }
[[ -d "$DB" ]]  || { echo "[mob_recon] DB missing: $DB"; exit 2; }

# Wenn Verzeichnis existiert → löschen, immer neu beginnen
if [[ -d "$OUT" ]]; then
  echo "[mob_recon] removing existing output dir $OUT"
  rm -rf "$OUT"
fi
mkdir -p "$OUT"

# Hauptbefehl
cmd=( mob_recon -i "$IN" -o "$OUT" -n "$CPUS" -d "$DB" --force )

echo "[mob_recon] ${cmd[*]}"
"${cmd[@]}"

date -u +"%Y-%m-%dT%H:%M:%SZ" > "$OUT/.done"
echo "[mob_recon] done -> $OUT"
