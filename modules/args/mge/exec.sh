#!/usr/bin/env bash
# modules/args/mge/exec.sh
# Mobile genetic element context for ARGs on all contigs (>= global.contig_min_len):
#   ISEScan          IS elements   -> mge/is_elements.tsv (contig, start, end, id, family)
#   IntegronFinder 2 integrons     -> mge/integrons.tsv   (contig, start, end, id, type)
# scripts/arg_context.py turns them into 'near_IS' (ARG within modules.args.is_window
# bp of an IS element) and 'in_integron'. Switch off with modules.args.mge: false.
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"
OUT="$ROOT/SRA/args/$SAMPLE/mge"; mkdir -p "$OUT"
DONE="$OUT/.done"
[[ "$FORCE" == "1" ]] && rm -f "$DONE"
[[ -s "$DONE" ]] && { echo "[mge] cached"; exit 0; }

# shellcheck source=bin/lib/cfg.sh
source "$ROOT/bin/lib/cfg.sh"
case "$(cfg modules.args.mge true)" in
  false|False|0|no) echo "[mge] disabled (modules.args.mge: false)"; date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"; exit 0 ;;
esac
MINLEN="$(cfg global.contig_min_len 1000)"
IN="$ROOT/SRA/assemblies/spades/$SAMPLE/contigs.len${MINLEN}.fasta"
[[ -s "$IN" ]] || { echo "[mge] no contigs >= ${MINLEN} bp: $IN"; exit 2; }

WORK="$OUT/work"; rm -rf "$WORK"; mkdir -p "$WORK"
cp -f "$IN" "$WORK/contigs.fasta"      # both tools name their outputs after the input file

# ---- IntegronFinder 2 ----
( cd "$WORK" && integron_finder --local-max --cpu "$CPUS" --outdir "$WORK/integronfinder" contigs.fasta ) \
  > "$OUT/integronfinder.log" 2>&1
INT="$(find "$WORK/integronfinder" -name '*.integrons' -print -quit)"
{
  printf 'contig\tstart\tend\tintegron_id\ttype\n'
  # one span per (replicon, integron): min(pos_beg)..max(pos_end); type = complete | In0 | CALIN
  [[ -n "$INT" ]] && awk -F'\t' '!/^#/ && $1 ~ /^integron_/ {
      k=$2 "\t" $1; if (!(k in s) || $4<s[k]) s[k]=$4; if ($5>e[k]) e[k]=$5; t[k]=$11 }
    END { for (k in s) { split(k, a, "\t"); print a[1] "\t" s[k] "\t" e[k] "\t" a[1] "|" a[2] "\t" t[k] } }' "$INT" | sort -k1,1 -k2,2n
} > "$OUT/integrons.tsv"
echo "[mge] integrons: $(( $(wc -l < "$OUT/integrons.tsv") - 1 ))"

# ---- ISEScan ----
( cd "$WORK" && isescan.py --seqfile contigs.fasta --output "$WORK/isescan" --nthread "$CPUS" ) \
  > "$OUT/isescan.log" 2>&1
IST="$(find "$WORK/isescan" -name '*.fasta.tsv' -print -quit)"
{
  printf 'contig\tstart\tend\tis_id\tfamily\n'
  [[ -n "$IST" ]] && awk -F'\t' 'NR>1 && NF>4 {print $1 "\t" $4 "\t" $5 "\t" $1 "|IS_" $4 "_" $5 "\t" $2}' "$IST" | sort -k1,1 -k2,2n
} > "$OUT/is_elements.tsv"
echo "[mge] IS elements: $(( $(wc -l < "$OUT/is_elements.tsv") - 1 ))"

# keep the tools' own result files, drop the large intermediates
mkdir -p "$OUT/raw"
[[ -n "$INT" ]] && cp -f "$INT" "$OUT/raw/"
[[ -n "$IST" ]] && cp -f "$IST" "$OUT/raw/isescan.tsv"
rm -rf "$WORK"
date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
echo "[mge] done -> $OUT"
