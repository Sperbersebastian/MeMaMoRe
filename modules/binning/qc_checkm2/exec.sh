#!/usr/bin/env bash
# modules/binning/qc_checkm2/exec.sh
# CheckM2 with robust DB bootstrap (v1.1.0 layout)
set -euo pipefail

# remove after testing 
#FORCE=1

: "${ROOT:?}"; : "${SAMPLE:?}"
PARAMS_YAML="${PARAMS_YAML:-$ROOT/config/params.yaml}"
FORCE="${FORCE:-0}"
CPUS="${CPUS:-}"

BINS="$ROOT/SRA/binning/magscot/$SAMPLE/bins"
OUT="$ROOT/SRA/binning/checkm2/$SAMPLE"
LOG="$ROOT/logs/binning_checkm2_${SAMPLE}.log"
mkdir -p "$OUT" "$ROOT/logs"

shopt -s nullglob
bins=( "$BINS"/*.fa* )
[[ ${#bins[@]} -gt 0 ]] || { echo "[checkm2] no bins"; exit 0; }

if [[ -z "$CPUS" ]]; then
  if command -v micromamba >/dev/null 2>&1; then
    CPUS="$(micromamba run -n env_binning yq -r '.global.threads // 4' "$PARAMS_YAML" 2>/dev/null || echo 4)"
  else CPUS=4; fi
fi

if [[ "$FORCE" != "1" ]]; then
  if compgen -G "$OUT/*" >/dev/null; then
    echo "[checkm2] existing results in $OUT; skip (use FORCE=1)"; exit 0
  fi
else
  rm -rf "$OUT"; mkdir -p "$OUT"
fi

# --- DB bootstrap compatible with 1.1.0
DBROOT="${CHECKM2_DB_DIR:-$ROOT/refdata/checkm2}"
mkdir -p "$DBROOT"
LOCK="$DBROOT/.download.lock"

# after download, the actual path is DBROOT/CheckM2_database
resolve_dbdir() {
  if [[ -d "$DBROOT/CheckM2_database" ]]; then
    echo "$DBROOT/CheckM2_database"
  else
    echo "$DBROOT"
  fi
}

db_ready() {
  local d; d="$(resolve_dbdir)"
  find "$d" -type f -name "*.dmnd" -print -quit | grep -q .
}

ensure_db() {
  if db_ready; then return 0; fi
  if mkdir "$LOCK" 2>/dev/null; then
    echo "[checkm2] downloading database into $DBROOT"
    find "$DBROOT" -mindepth 1 -maxdepth 1 -exec rm -rf {} + 2>/dev/null || true
    micromamba run -n env_binning checkm2 database --download --path "$DBROOT"
    rmdir "$LOCK" || true
  else
    echo "[checkm2] waiting for DB download to finish..."
    for _ in {1..180}; do [[ -d "$LOCK" ]] || break; sleep 10; done
  fi
  db_ready || { echo "[checkm2] DB not ready in $DBROOT" >&2; exit 1; }
}

ensure_db
# resolve DB dir and dmnd file
DBDIR="$(resolve_dbdir)"                              # e.g. /.../refdata/checkm2/CheckM2_database
DBFILE="$(find "$DBDIR" -type f -name '*.dmnd' | head -n1)" || true
[[ -s "$DBFILE" ]] || { echo "[checkm2] no .dmnd in $DBDIR"; exit 1; }

export CHECKM2DB="$DBDIR"   # models live here; CheckM2 reads it internally

echo "[checkm2] start (threads=$CPUS) -> $OUT ; db=$DBFILE"
set +e
micromamba run -n env_binning checkm2 predict \
  --input "$BINS" \
  --output-directory "$OUT" \
  --database_path "$DBFILE" \
  -x fa \
  --threads "$CPUS" \
  --force \
  >"$LOG" 2>&1
rc=$?
set -e


if [[ $rc -ne 0 ]]; then
  echo "[checkm2] failed, see $LOG"; exit $rc
fi
echo "[checkm2] -> $OUT"
