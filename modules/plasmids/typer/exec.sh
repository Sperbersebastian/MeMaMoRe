#!/usr/bin/env bash
# MOB-typer + optional geNomad on dereplicated plasmids (multi mode)
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

IN="$ROOT/SRA/plasmids/$SAMPLE/cluster/plasmids_derep.$SAMPLE.fasta"
OUT="$ROOT/SRA/plasmids/$SAMPLE/typer"; mkdir -p "$OUT"

[[ -s "$IN" ]] || { echo "[mob_typer] no input: $IN"; exit 1; }
command -v mob_typer >/dev/null 2>&1 || { echo "[mob_typer] mob_typer not in PATH"; exit 2; }

# DBs if present
MOBS_DB_DEFAULT="$ROOT/refdata/mobsuite_db"
GENOMAD_DB_DEFAULT="$ROOT/refdata/genomad/genomad_db"
[[ -d "$MOBS_DB_DEFAULT" ]]    && MOBS_DB="$MOBS_DB_DEFAULT"       || MOBS_DB=""
[[ -d "$GENOMAD_DB_DEFAULT" ]] && GENOMAD_DB="$GENOMAD_DB_DEFAULT" || GENOMAD_DB=""

# cache flag
DONE="$OUT/.done"
[[ "$FORCE" == "1" ]] && rm -f "$DONE"
[[ -s "$DONE" ]] && { echo "[mob_typer] cached"; exit 0; }

# --- MOB-typer (v3.1.9) in multi mode ---
OUTFILE="$OUT/mobtyper_results.txt"
CMD=(mob_typer -i "$IN" -o "$OUTFILE" -n "$CPUS" --multi)
[[ -n "$MOBS_DB" ]] && CMD+=(-d "$MOBS_DB")
[[ -n "${MOB_TYPER_ARGS:-}" ]] && CMD+=(${MOB_TYPER_ARGS})

echo "[mob_typer] ${CMD[*]}"
"${CMD[@]}"

# stable symlink
[[ -s "$OUTFILE" && ! -e "$OUT/results.mob_typer.tsv" ]] && ln -sf "$(basename "$OUTFILE")" "$OUT/results.mob_typer.tsv" || true

# --- geNomad end-to-end on clustered FASTA (default on) ---
GENOMAD_ON_CLUSTERED="${GENOMAD_ON_CLUSTERED:-1}"
if [[ "$GENOMAD_ON_CLUSTERED" == "1" ]]; then
  if [[ -n "$GENOMAD_DB" ]]; then
    GOUT="$OUT/genomad_on_derep"; mkdir -p "$GOUT"
    echo "[genomad] end-to-end on dereplicated plasmids"
    micromamba run -n genomad_env genomad end-to-end "$IN" "$GOUT" "$GENOMAD_DB" --threads "$CPUS"
  else
    echo "[genomad] DB missing at $GENOMAD_DB_DEFAULT -> skip"
  fi
fi

date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
echo "[mob_typer] done -> $OUT"
