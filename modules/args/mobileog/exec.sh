#!/usr/bin/env bash
# modules/args/mobileog/exec.sh
# mobileOG-db (mobile genetic element protein families) via DIAMOND blastp on
# the central gene calls (genes.faa). Output SRA/args/<s>/mobileog/genes.tsv:
# gene_id, mobileOG id, name, major category (phage, integration/excision,
# replication/recombination/repair, transfer, stability/transfer/defense),
# minor categories, pident, qcov, evalue. scripts/arg_context.py reports the
# mobileOG categories of genes within modules.args.is_window of each ARG.
# DB: refdata/mobileog/mobileOGdb_90-2.0.fasta (mobileOG-db 2.0.1, Zenodo
# 14725479, CC BY 4.0, 2.3 GB), DIAMOND index built once next to it.
# Thresholds as in mobileOGs-pl: identity 90, query cover 90, e-value 1e-20.
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"
OUT="$ROOT/SRA/args/$SAMPLE/mobileog"; mkdir -p "$OUT"
DONE="$OUT/.done"
[[ "$FORCE" == "1" ]] && rm -f "$DONE"
[[ -s "$DONE" ]] && { echo "[mobileog] cached"; exit 0; }

# shellcheck source=bin/lib/cfg.sh
source "$ROOT/bin/lib/cfg.sh"
case "$(cfg modules.args.mobileog true)" in
  false|False|0|no) echo "[mobileog] disabled (modules.args.mobileog: false)"; date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"; exit 0 ;;
esac
ID="$(cfg modules.args.mobileog_min_id 90)"
QCOV="$(cfg modules.args.mobileog_min_qcov 90)"
EVALUE="$(cfg modules.args.mobileog_evalue 1e-20)"

GENES="$ROOT/SRA/assemblies/genes/$SAMPLE/genes.faa"
[[ -s "$GENES" ]] || { echo "[mobileog] missing $GENES (run assembly)"; exit 2; }
command -v diamond >/dev/null 2>&1 || { echo "[mobileog] diamond not in PATH"; exit 3; }

DBDIR="${MOBILEOG_DIR:-$ROOT/refdata/mobileog}"
FASTA="$DBDIR/mobileOGdb_90-2.0.fasta"
DMND="$DBDIR/mobileOGdb_90-2.0.dmnd"
if [[ ! -s "$FASTA" ]]; then
  echo "[mobileog] DB missing -> downloading mobileOG-db 2.0.1 (2.3 GB) into $DBDIR"
  mkdir -p "$DBDIR"
  curl -fL -C - -o "$FASTA.part" "https://zenodo.org/records/14725479/files/mobileOGdb_90-2.0.fasta?download=1"
  mv -f "$FASTA.part" "$FASTA"
fi
if [[ ! -s "$DMND" ]]; then
  echo "[mobileog] building DIAMOND index $DMND"
  diamond makedb --in "$FASTA" -d "${DMND%.dmnd}" -p "$CPUS" --quiet
fi

sed '/^>/!s/\*//g' "$GENES" > "$OUT/genes.faa"
diamond blastp -q "$OUT/genes.faa" -d "$DMND" -o "$OUT/hits.m8" -p "$CPUS" --quiet \
  -k 1 --id "$ID" --query-cover "$QCOV" -e "$EVALUE" \
  --outfmt 6 qseqid sseqid pident qcovhsp evalue bitscore
rm -f "$OUT/genes.faa"

# subject IDs: mobileOG_id|source_id|name|major category|minor categories|source|evidence
{
  printf 'gene_id\tmobileog_id\tname\tmajor_category\tminor_categories\tpident\tqcov\tevalue\n'
  awk -F'\t' 'BEGIN{OFS="\t"} { n=split($2, h, "|"); print $1, h[1], (n>2?h[3]:""), (n>3?h[4]:""), (n>4?h[5]:""), $3, $4, $5 }' "$OUT/hits.m8"
} > "$OUT/genes.tsv"
rm -f "$OUT/hits.m8"
echo "[mobileog] $(( $(wc -l < "$OUT/genes.tsv") - 1 )) genes with a mobileOG hit -> $OUT/genes.tsv"
date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
