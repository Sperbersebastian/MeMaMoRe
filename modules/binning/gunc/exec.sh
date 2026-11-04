#!/usr/bin/env bash
# modules/binning/gunc/exec.sh
# GUNC chimerism check on the refined MAGScoT bins, using the central gene calls
# (genes.faa split per bin, --gene_calls). Result: SRA/binning/gunc/<sample>/
# GUNC.<db>.maxCSS_level.tsv (column pass.GUNC), read by scripts/summarize_samples.py.
# env in: ROOT SAMPLE [CPUS] [FORCE] [PARAMS_YAML] [GUNC_DB]
set -euo pipefail

: "${ROOT:?}"; : "${SAMPLE:?}"
CPUS="${CPUS:-8}"; FORCE="${FORCE:-0}"
OUT="$ROOT/SRA/binning/gunc/$SAMPLE"
BINS="$ROOT/SRA/binning/magscot/$SAMPLE/bins"
MAP="$ROOT/SRA/binning/magscot/$SAMPLE/contigs_to_bin.tsv"
GENES="$ROOT/SRA/assemblies/genes/$SAMPLE/genes.faa"
LOG="$ROOT/logs/binning_gunc_${SAMPLE}.log"
mkdir -p "$OUT" "$ROOT/logs"

cfg(){ local v=""; [[ -n "${PARAMS_YAML:-}" && -s "$PARAMS_YAML" ]] && \
         v="$(micromamba run -n env_sra_tools python "$ROOT/scripts/merge_params.py" --defaults "$PARAMS_YAML" --get "$1" 2>/dev/null || true)"
       [[ -n "$v" && "$v" != "None" && "$v" != "null" ]] && echo "$v" || echo "$2"; }
case "$(cfg modules.binning.gunc true)" in
  false|False|0|no) echo "[gunc] disabled (modules.binning.gunc: false)"; exit 0 ;;
esac

shopt -s nullglob; bins=( "$BINS"/*.fa ); shopt -u nullglob
[[ ${#bins[@]} -gt 0 ]] || { echo "[gunc] no refined bins for $SAMPLE -> skip"; exit 0; }
[[ -s "$GENES" && -s "$MAP" ]] || { echo "[gunc] missing $GENES or $MAP -> skip"; exit 0; }

RES="$(find "$OUT" -maxdepth 1 -name 'GUNC.*.maxCSS_level.tsv' -print -quit 2>/dev/null || true)"
if [[ -n "$RES" && "$FORCE" != "1" && "$RES" -nt "$MAP" ]]; then
  echo "[gunc] cached -> $RES"; exit 0
fi

# env + reference DB (ProGenomes 2.1 by default, ~13 GB, fetched once)
micromamba env list --json | grep -q '/gunc_env"' || bash "$ROOT/envs/create_env_gunc.sh"
DBDIR="$ROOT/refdata/gunc"
DB="${GUNC_DB:-$(find "$DBDIR" -maxdepth 2 -name 'gunc_db_*.dmnd' 2>/dev/null | sort | head -n1 || true)}"
if [[ -z "$DB" || ! -s "$DB" ]]; then
  echo "[gunc] DB missing -> downloading ProGenomes 2.1 into $DBDIR"
  mkdir -p "$DBDIR"
  micromamba run -n gunc_env gunc download_db "$DBDIR" -db progenomes_2.1 >>"$LOG" 2>&1
  DB="$(find "$DBDIR" -maxdepth 2 -name 'gunc_db_*.dmnd' | sort | head -n1)"
fi
[[ -s "$DB" ]] || { echo "[gunc] ERROR: no GUNC DB in $DBDIR"; exit 2; }

# proteins per bin from the central gene calls (gene IDs are <contig>_<n>)
FAA_DIR="$OUT/faa"; rm -rf "$FAA_DIR" "$OUT"/GUNC.* "$OUT/gunc_output" "$OUT/diamond_output"; mkdir -p "$FAA_DIR"
awk -v dir="$FAA_DIR" '
  FNR==NR { if (FNR>1) bin[$1]=$2; next }
  /^>/ { id=substr($1,2); c=id; sub(/_[0-9]+$/, "", c); out=(c in bin) ? dir "/" bin[c] ".faa" : "" }
  out!="" { print > out }
' "$MAP" "$GENES"

echo "[gunc] ${#bins[@]} bins, db=$(basename "$DB")"
micromamba run -n gunc_env gunc run -d "$FAA_DIR" -e .faa -g -r "$DB" -t "$CPUS" -o "$OUT" >>"$LOG" 2>&1
RES="$(find "$OUT" -maxdepth 1 -name 'GUNC.*.maxCSS_level.tsv' -print -quit)"
[[ -s "$RES" ]] || { echo "[gunc] ERROR: no GUNC result in $OUT (see $LOG)"; exit 3; }
echo "[gunc] $(awk -F'\t' 'NR>1 && $NF=="True"' "$RES" | wc -l)/$(( $(wc -l < "$RES") - 1 )) bins pass -> $RES"
