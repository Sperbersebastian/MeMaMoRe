#!/usr/bin/env bash
set -euo pipefail
# env in: ROOT SAMPLE PARAMS_YAML [FORCE] [CPUS]
BINS="$ROOT/SRA/binning/magscot/$SAMPLE/bins"
OUT="$ROOT/SRA/binning/checkm2/$SAMPLE"
LOG="$ROOT/logs/binning_checkm2_${SAMPLE}.log"
mkdir -p "$OUT" "$ROOT/logs"

# bins present?
shopt -s nullglob
bins=("$BINS"/*.fa*); [[ ${#bins[@]} -gt 0 ]] || { echo "[checkm2] no bins"; exit 0; }

# threads (fallback to 4)
get_threads(){ micromamba run -n env_sra_tools python "$ROOT/scripts/merge_params.py" \
  --defaults "$PARAMS_YAML" --get global.threads 2>/dev/null || echo 4; }
CPUS="${CPUS:-$(get_threads)}"

# force/skip logic
if [[ "${FORCE:-0}" != "1" ]]; then
  if compgen -G "$OUT/*" >/dev/null; then
    echo "[checkm2] existing results in $OUT; skip (use --force to overwrite)"
    exit 0
  fi
else
  rm -rf "$OUT"
  mkdir -p "$OUT"
fi

# run
micromamba run -n env_binning checkm2 predict \
  --genome_dir "$BINS" -x fa \
  --out_dir "$OUT" \
  --threads "$CPUS" \
  --force \
  >"$LOG" 2>&1 || true

echo "[checkm2] -> $OUT"
