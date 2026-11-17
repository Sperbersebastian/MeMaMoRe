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
[[ -d "$BINS" ]] || { echo "[gtdbtk] no refined bins ($BINS missing) -> skip"; exit 0; }
shopt -s nullglob
fa_files=( "$BINS"/*.fa )
echo "[gtdbtk] .fa files found: ${#fa_files[@]} in $BINS"
[[ ${#fa_files[@]} -gt 0 ]] || { echo "[gtdbtk] no .fa bins -> exit 0"; exit 0; }

# --- memory preflight: classify_wf needs ~64 GB (pplacer; RAM + swap count,
#     pplacer also uses --scratch_dir); below that skip with a clear message
# shellcheck source=bin/lib/memory.sh
source "$ROOT/bin/lib/memory.sh"
mem_gate gtdbtk "${GTDBTK_MIN_MEM_GB:-64}" || exit 0

# --- env: GTDB-Tk lives in its own env (gtdbtk_env, envs/create_env_gtdbtk.sh)
micromamba env list --json | grep -q '/gtdbtk_env"' || bash "$ROOT/envs/create_env_gtdbtk.sh"
GTDB_ENV="gtdbtk_env"

# --- DB: newest refdata/gtdbtk/release<N> (or GTDBTK_DATA_PATH), verified via check_install
DBROOT_INPUT="${GTDBTK_DATA_PATH:-$ROOT/refdata/gtdbtk}"
resolve_dbroot() {
  local d="$1"
  [[ -d "$d/markers" ]] && { echo "$d"; return; }
  local cand
  cand="$(find "$d" -maxdepth 2 -type d -name 'release*' -exec test -d '{}/markers' \; -print | sort -V | tail -n1 || true)"
  [[ -n "${cand:-}" ]] && { echo "$cand"; return; }
  echo "$d"
}
DBROOT="$(resolve_dbroot "$DBROOT_INPUT")"
export GTDBTK_DATA_PATH="$DBROOT"
echo "[gtdbtk] DB candidate: $GTDBTK_DATA_PATH"

if ! micromamba run -n "$GTDB_ENV" gtdbtk check_install >/dev/null; then
  echo "[gtdbtk] ERROR: gtdbtk check_install failed for: $GTDBTK_DATA_PATH"
  echo "[gtdbtk] Hints:"
  echo "  - the path must be the release folder holding markers/, pplacer/, taxonomy/ ..."
  echo "  - GTDB-Tk >= 2.7 needs release232 or newer: scripts/setup_gtdbtk_data.sh latest"
  echo "  - content:"
  ls -1 "$GTDBTK_DATA_PATH" || true
  exit 1
fi
echo "[gtdbtk] DB verified by check_install."

# --- skip if results exist unless FORCE=1
# refined bins identity (magscot/<s>/bins/.md5) + reference release: results are reused
# only for the same bins classified against the same GTDB release
BINS_SUM="$(cat "$ROOT/SRA/binning/magscot/$SAMPLE/bins/.md5" 2>/dev/null || echo none)"
STAMP="$BINS_SUM.$(basename "$DBROOT")"
if [[ "$FORCE" != "1" && "$(cat "$OUT/.bins.md5" 2>/dev/null)" == "$STAMP" ]] && { [[ -s "$OUT/gtdbtk.bac120.summary.tsv" ]] || [[ -s "$OUT/gtdbtk.ar53.summary.tsv" ]]; }; then
  echo "[gtdbtk] outputs exist and FORCE!=1 -> skip"
  echo "[gtdbtk] -> $OUT"
  echo "[gtdbtk] log -> $LOGFILE"
  exit 0
fi

# --- stale state: GTDB-Tk resumes from $OUT/gtdbtk.json and crashes on a state file of
#     another version (2.8: "unexpected keyword argument 'skip_ani_screen'"); intermediates
#     of other bins or another release must not be reused either. Keep them only for an
#     interrupted run of the same bins, release and GTDB-Tk version.
GTDB_VER="$(micromamba run -n "$GTDB_ENV" gtdbtk --version 2>/dev/null | awk '/version/{print $3; exit}' || true)"
RUN_ID="$STAMP.${GTDB_VER:-unknown}"
if [[ "$FORCE" == "1" || "$(cat "$OUT/.run.id" 2>/dev/null)" != "$RUN_ID" ]]; then
  echo "[gtdbtk] clearing previous GTDB-Tk output in $OUT"
  rm -rf "$OUT/identify" "$OUT/align" "$OUT/classify" "$OUT/gtdbtk.json" "$OUT/gtdbtk.log" \
         "$OUT/gtdbtk.warnings.log" "$OUT/.bins.md5"
  find "$OUT" -maxdepth 1 -name 'gtdbtk.*.tsv' -delete
fi
echo "$RUN_ID" > "$OUT/.run.id"

# --- run classify_wf on the .fa bins
# --skip_ani_screen exists only up to GTDB-Tk 2.6 (later versions always screen with skani)
ANI_OPT=()
micromamba run -n "$GTDB_ENV" gtdbtk classify_wf -h 2>&1 | grep -q -- '--skip_ani_screen' && ANI_OPT=(--skip_ani_screen)
echo "[gtdbtk] RUN classify_wf cpus=$CPUS bins_dir=$BINS extension=fa"
micromamba run -n "$GTDB_ENV" gtdbtk classify_wf \
  --genome_dir "$BINS" \
  --out_dir "$OUT" \
  --extension fa \
  --cpus "$CPUS" \
  --pplacer_cpus "$CPUS" \
  --scratch_dir "$SCRATCH" \
  --tmpdir "$TMPDIR" \
  "${ANI_OPT[@]}"

echo "$STAMP" > "$OUT/.bins.md5"
echo "[gtdbtk] DONE -> $OUT"
echo "[gtdbtk] LOG  -> $LOGFILE"
