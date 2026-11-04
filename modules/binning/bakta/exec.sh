#!/usr/bin/env bash
# modules/binning/bakta/exec.sh
# Optional Bakta annotation of the refined MAGs (off by default:
# modules.binning.bakta: true to enable). DB type modules.binning.bakta_db:
# light (default, ~1.4 GB packed) or full (~37 GB packed), fetched once into
# refdata/bakta/<type>. Output SRA/binning/bakta/<s>/<bin>/<bin>.{gff3,gbff,faa,tsv,...};
# reused while the refined bins are unchanged (magscot bins/.md5).
# env in: ROOT SAMPLE [CPUS] [FORCE] [PARAMS_YAML]
set -euo pipefail

: "${ROOT:?}"; : "${SAMPLE:?}"
CPUS="${CPUS:-8}"; FORCE="${FORCE:-0}"
# shellcheck source=bin/lib/cfg.sh
source "$ROOT/bin/lib/cfg.sh"
case "$(cfg modules.binning.bakta false)" in
  true|True|1|yes) ;;
  *) exit 0 ;;                      # default: off, silent
esac
DBTYPE="$(cfg modules.binning.bakta_db light)"
[[ "$DBTYPE" == light || "$DBTYPE" == full ]] || { echo "[bakta] invalid modules.binning.bakta_db '$DBTYPE' (light|full)"; exit 2; }

BINS="$ROOT/SRA/binning/magscot/$SAMPLE/bins"
OUT="$ROOT/SRA/binning/bakta/$SAMPLE"
LOG="$ROOT/logs/binning_bakta_${SAMPLE}.log"
mkdir -p "$OUT" "$ROOT/logs"
shopt -s nullglob; bins=( "$BINS"/*.fa ); shopt -u nullglob
[[ ${#bins[@]} -gt 0 ]] || { echo "[bakta] no refined bins for $SAMPLE -> skip"; exit 0; }

BINS_SUM="$(cat "$BINS/.md5" 2>/dev/null || echo none)"
if [[ "$FORCE" != "1" && "$(cat "$OUT/.bins.md5" 2>/dev/null)" == "$BINS_SUM.$DBTYPE" ]]; then
  echo "[bakta] results for these bins exist in $OUT; skip"; exit 0
fi

micromamba env list --json | grep -q '/bakta_env"' || bash "$ROOT/envs/create_env_bakta.sh"
DBROOT="${BAKTA_DB_DIR:-$ROOT/refdata/bakta/$DBTYPE}"
DB="$(find "$DBROOT" -maxdepth 2 -name version.json -printf '%h\n' 2>/dev/null | head -n1 || true)"
if [[ -z "$DB" ]]; then
  echo "[bakta] DB ($DBTYPE) missing -> downloading into $DBROOT"
  mkdir -p "$DBROOT"
  micromamba run -n bakta_env bakta_db download --output "$DBROOT" --type "$DBTYPE" >>"$LOG" 2>&1
  DB="$(find "$DBROOT" -maxdepth 2 -name version.json -printf '%h\n' | head -n1)"
fi
[[ -n "$DB" ]] || { echo "[bakta] ERROR: no Bakta DB in $DBROOT"; exit 3; }

rm -rf "$OUT"; mkdir -p "$OUT"
for b in "${bins[@]}"; do
  n="$(basename "$b" .fa)"
  echo "[bakta] $SAMPLE/$n" | tee -a "$LOG"
  micromamba run -n bakta_env bakta --db "$DB" --meta --threads "$CPUS" --skip-plot --force \
    --prefix "$n" --output "$OUT/$n" "$b" >>"$LOG" 2>&1
done
echo "$BINS_SUM.$DBTYPE" > "$OUT/.bins.md5"
echo "[bakta] ${#bins[@]} MAG(s) annotated ($DBTYPE DB) -> $OUT"
