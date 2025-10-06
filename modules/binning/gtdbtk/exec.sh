#!/usr/bin/env bash
set -euxo pipefail
# GTDB-Tk classify wrapper
# Required: ROOT, SAMPLE
# Optional: CPUS, FORCE, GTDBTK_DATA_PATH, SCRATCH_DIR, TMPDIR

: "${ROOT:?ROOT not set}"
: "${SAMPLE:?SAMPLE not set}"

# --- logging (line-buffered)
LOGDIR="$ROOT/logs"; mkdir -p "$LOGDIR"
LOGFILE="$LOGDIR/binning_gtdb_tk_${SAMPLE}.log"
exec > >(stdbuf -oL tee -a "$LOGFILE") 2>&1
echo "[gtdbtk] START sample=$SAMPLE root=$ROOT host=$(hostname) user=$(whoami)"

# --- paths
OUT="$ROOT/SRA/binning/gtdbtk/$SAMPLE"; mkdir -p "$OUT"
BINS="$ROOT/SRA/binning/magscot/$SAMPLE/bins"
CPUS="${CPUS:-8}"
FORCE="${FORCE:-0}"
SCRATCH="${SCRATCH_DIR:-$OUT/.scratch}"
TMPDIR_IN="${TMPDIR:-$OUT/.tmp}"
mkdir -p "$SCRATCH" "$TMPDIR_IN"
export TMPDIR="$TMPDIR_IN"

# --- inputs
[[ -d "$BINS" ]] || { echo "[gtdbtk] ERROR: bins dir not found: $BINS"; exit 1; }
shopt -s nullglob
fa_files=( "$BINS"/*.fa )
echo "[gtdbtk] .fa files found: ${#fa_files[@]} in $BINS"
[[ ${#fa_files[@]} -gt 0 ]] || { echo "[gtdbtk] no .fa bins -> exit 0"; exit 0; }

# --- DB: resolve + verify via check_install
DBROOT_INPUT="${GTDBTK_DATA_PATH:-$ROOT/refdata/gtdbtk/release226}"
resolve_dbroot() {
  local d="$1"
  [[ -d "$d/markers" ]] && { echo "$d"; return; }
  [[ -d "$d/release226/markers" ]] && { echo "$d/release226"; return; }
  local cand
  cand="$(find "$d" -maxdepth 2 -type d -name 'release*' -exec test -d '{}/markers' \; -print | head -n1 || true)"
  [[ -n "${cand:-}" ]] && { echo "$cand"; return; }
  echo "$d"
}
DBROOT="$(resolve_dbroot "$DBROOT_INPUT")"
export GTDBTK_DATA_PATH="$DBROOT"
echo "[gtdbtk] DB candidate: $GTDBTK_DATA_PATH"

# Verifiziere DB mit GTDB-Tk selbst
if ! micromamba run -n env_binning gtdbtk check_install >/dev/null; then
  echo "[gtdbtk] ERROR: gtdbtk check_install failed for: $GTDBTK_DATA_PATH"
  echo "[gtdbtk] Hints:"
  echo "  - Pfad muss auf den *Release*-Ordner zeigen, der 'markers/', 'pplacer/', 'taxonomy/' etc. enthält."
  echo "  - Beispiel: /media/Box/MeMaMoRe/refdata/gtdbtk/release226/"
  echo "  - Inhalt:"
  ls -1 "$GTDBTK_DATA_PATH" || true
  exit 1
fi
echo "[gtdbtk] DB verified by check_install."

# --- skip if results exist unless FORCE=1
if [[ "$FORCE" != "1" ]] && { [[ -s "$OUT/gtdbtk.bac120.summary.tsv" ]] || [[ -s "$OUT/gtdbtk.ar53.summary.tsv" ]]; }; then
  echo "[gtdbtk] outputs exist and FORCE!=1 -> skip"
  echo "[gtdbtk] -> $OUT"
  echo "[gtdbtk] log -> $LOGFILE"
  exit 0
fi

# --- run classify_wf auf .fa-Bins
echo "[gtdbtk] RUN classify_wf cpus=$CPUS bins_dir=$BINS extension=fa"
micromamba run -n env_binning gtdbtk classify_wf \
  --genome_dir "$BINS" \
  --out_dir "$OUT" \
  --extension fa \
  --cpus "$CPUS" \
  --pplacer_cpus "$CPUS" \
  --scratch_dir "$SCRATCH" \
  --tmpdir "$TMPDIR" \
  --skip_ani_screen

echo "[gtdbtk] DONE -> $OUT"
echo "[gtdbtk] LOG  -> $LOGFILE"
