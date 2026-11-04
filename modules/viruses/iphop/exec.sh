#!/usr/bin/env bash
# modules/viruses/iphop/exec.sh
# iPHoP host prediction for the virus contigs that CheckV rates at least
# medium-quality (all union contigs if CheckV is missing).
# Optional step: bin/modules/viruses.sh only calls it when the iPHoP DB exists
# (~340 GB, `iphop download --db_dir "$ROOT/refdata/iphop"`).
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
# Env:   IPHOP_DB (default: first $ROOT/refdata/iphop/*_pub_rw dir)
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

OUT="$ROOT/SRA/viruses/$SAMPLE/iphop"; mkdir -p "$OUT"
DONE="$OUT/.done"
[[ "$FORCE" == "1" ]] && rm -f "$DONE"
[[ -s "$DONE" ]] && { echo "[iphop] cached"; exit 0; }

command -v iphop >/dev/null 2>&1 || { echo "[iphop] iphop not in PATH"; exit 3; }

DB="${IPHOP_DB:-$(find "$ROOT/refdata/iphop" -mindepth 1 -maxdepth 1 -type d -name '*_rw' 2>/dev/null | sort | tail -n1)}"
[[ -n "$DB" && -d "$DB" ]] || { echo "[iphop] DB missing (set IPHOP_DB or download into $ROOT/refdata/iphop)"; exit 2; }

UNION="$ROOT/SRA/viruses/$SAMPLE/union/viruses_union.fasta"
CHECKV="$ROOT/SRA/viruses/$SAMPLE/checkv/quality_summary.tsv"
[[ -s "$UNION" ]] || { echo "[iphop] no input: $UNION"; exit 1; }

# keep Complete / High-quality / Medium-quality contigs (iPHoP is slow on fragments)
IN="$OUT/input_viruses.fasta"
if [[ -s "$CHECKV" ]]; then
  awk -F'\t' 'NR==FNR { if (FNR>1 && $8 ~ /^(Complete|High-quality|Medium-quality)$/) keep[$1]=1; next }
              /^>/ { id=substr($1,2); p=(id in keep) } p' "$CHECKV" "$UNION" > "$IN"
else
  cp -f "$UNION" "$IN"
fi
n=$(grep -c '^>' "$IN" || true)
[[ "$n" -gt 0 ]] || { echo "[iphop] no medium+ quality viruses for $SAMPLE"; date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"; exit 0; }

echo "[iphop] predicting hosts for $n virus contig(s) with $DB"
iphop predict --fa_file "$IN" --db_dir "$DB" --out_dir "$OUT/result" --num_threads "$CPUS"

# genus-level table (score >= 90) for the GUI and the multi-sample summary
G="$(find "$OUT/result" -maxdepth 1 -name 'Host_prediction_to_genus_m*.csv' | head -n1 || true)"
if [[ -n "$G" ]]; then
  ln -sf "${G#"$OUT"/}" "$OUT/host_genus.csv"
  echo "[iphop] $(($(wc -l < "$G") - 1)) genus-level prediction(s) -> $OUT/host_genus.csv"
else
  echo "[iphop] WARNING: no genus-level prediction file in $OUT/result"
fi
date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
echo "[iphop] done -> $OUT"
