#!/usr/bin/env bash
# modules/args/summary/exec.sh
# Merge ARG results from DeepARG, RGI, AMR++ across all genomic contexts
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

OUT="$ROOT/SRA/args/$SAMPLE/summary"; mkdir -p "$OUT"
DONE="$OUT/.done"
[[ "$FORCE" == "1" ]] && rm -f "$DONE"
[[ -s "$DONE" ]] && { echo "[arg-summary] cached"; exit 0; }

SUMMARY="$OUT/args_summary.tsv"

# Header
echo -e "source\tcontext\tgene\tARG_class\tcontig_or_gene\tidentity\tevalue_or_score\tnotes" > "$SUMMARY"

count=0

# --- DeepARG (all contexts) ---
DEEPARG_DIR="$ROOT/SRA/args/$SAMPLE/deeparg"
for ctx_dir in "$DEEPARG_DIR"/contigs "$DEEPARG_DIR"/mags "$DEEPARG_DIR"/plasmids "$DEEPARG_DIR"/viruses; do
  [[ -d "$ctx_dir" ]] || continue
  context=$(basename "$ctx_dir")
  for f in "$ctx_dir"/*.mapping.ARG; do
    [[ -s "$f" ]] || continue
    label=$(basename "$f" .mapping.ARG)
    echo "[arg-summary] adding DeepARG/$context/$label"
    tail -n +2 "$f" | awk -F'\t' -v ctx="$context" -v lbl="$label" \
      'BEGIN{OFS="\t"} {print "deeparg",ctx,$1,$6,$2,$4,$5,"prob="$8";label="lbl}' >> "$SUMMARY"
    ((count++)) || true
  done
done

# --- RGI (all contexts) ---
RGI_DIR="$ROOT/SRA/args/$SAMPLE/rgi"
for ctx_dir in "$RGI_DIR"/contigs "$RGI_DIR"/mags "$RGI_DIR"/plasmids "$RGI_DIR"/viruses; do
  [[ -d "$ctx_dir" ]] || continue
  context=$(basename "$ctx_dir")
  for f in "$ctx_dir"/*.txt; do
    [[ -s "$f" ]] || continue
    label=$(basename "$f" .txt)
    echo "[arg-summary] adding RGI/$context/$label"
    tail -n +2 "$f" | awk -F'\t' -v ctx="$context" -v lbl="$label" \
      'BEGIN{OFS="\t"} {print "rgi",ctx,$9,$15,$1,$10,$11,"model="$12";label="lbl}' >> "$SUMMARY"
    ((count++)) || true
  done
done

# --- AMR++ (reads only) ---
AMRPP_FILE="$ROOT/SRA/args/$SAMPLE/amrplusplus/resistome_counts.tsv"
if [[ -s "$AMRPP_FILE" ]]; then
  echo "[arg-summary] adding AMR++/reads"
  tail -n +2 "$AMRPP_FILE" | awk -F'\t' \
    'BEGIN{OFS="\t"} {print "amrplusplus","reads",$1,"megares",$1,"NA","NA","reads="$3}' >> "$SUMMARY"
  ((count++)) || true
fi

[[ $count -gt 0 ]] || echo "[arg-summary] WARNING: no ARG results from any tool"

TOTAL=$(tail -n +2 "$SUMMARY" | wc -l)
echo "[arg-summary] unified table: $TOTAL entries from $count sources"

# --- Per-context summary ---
echo "[arg-summary] breakdown by context:"
tail -n +2 "$SUMMARY" | awk -F'\t' '{ctx[$2]++; tool[$1]++} END{
  for(c in ctx) printf "  %s: %d entries\n",c,ctx[c]
  for(t in tool) printf "  [%s]: %d entries\n",t,tool[t]
}'

date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
echo "[arg-summary] done -> $SUMMARY"
