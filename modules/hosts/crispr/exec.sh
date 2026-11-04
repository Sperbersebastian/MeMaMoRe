#!/usr/bin/env bash
# modules/hosts/crispr/exec.sh
# CRISPR spacers of the refined MAGs (MinCED) matched against the plasmid and
# virus contigs of all samples (blastn-short): a spacer hit is host evidence
# (the MAG's CRISPR system has seen that element).
# Targets: global plasmid cluster representatives ("<sample>|<tool>|<contig>")
#          and every sample's virus union (prefixed "<sample>|").
# A hit needs full spacer length (qcov 100) and <= modules.hosts.crispr_max_mismatch
# mismatches + gaps (default 1). Output SRA/summary/crispr_hits.tsv
# (mag, spacer, target, target_type, mismatches, identity) + SRA/catalog/crispr/.
# Runs inside crispr_env (MinCED, BLAST). env in: ROOT SAMPLES (space-separated) [PARAMS_YAML]
set -euo pipefail
: "${ROOT:?}"; : "${SAMPLES:?}"
# shellcheck source=bin/lib/cfg.sh
source "$ROOT/bin/lib/cfg.sh"
CPUS="$(cfg global.threads 8)"
MAXMM="$(cfg modules.hosts.crispr_max_mismatch 1)"
OUT="$ROOT/SRA/catalog/crispr"; rm -rf "$OUT"; mkdir -p "$OUT/minced"
TSV="$ROOT/SRA/summary/crispr_hits.tsv"
mkdir -p "$(dirname "$TSV")"
read -r -a S <<<"$SAMPLES"

# ---- spacers of every refined MAG ("<sample>__<bin>|spacerN") ----
SP="$OUT/spacers.fa"; : > "$SP"
for s in "${S[@]}"; do
  shopt -s nullglob
  for b in "$ROOT/SRA/binning/magscot/$s/bins"/*.fa; do
    m="${s}__$(basename "$b" .fa)"
    minced -spacers "$b" "$OUT/minced/$m.crisprs" "$OUT/minced/$m.gff" >/dev/null 2>&1 || continue
    [[ -s "$OUT/minced/${m}_spacers.fa" ]] && awk -v m="$m" '/^>/{print ">" m "|" substr($1,2); next} {print}' \
      "$OUT/minced/${m}_spacers.fa" >> "$SP"
  done
  shopt -u nullglob
done
NSP=$(grep -c '^>' "$SP" || true)

# ---- targets ----
TG="$OUT/targets.fa"; : > "$TG"
GPL="$ROOT/SRA/plasmids/_global/cluster/plasmids_derep.global.fasta"
[[ -s "$GPL" ]] && awk '/^>/{print ">plasmid|" substr($1,2); next} {print}' "$GPL" >> "$TG"
for s in "${S[@]}"; do
  V="$ROOT/SRA/viruses/$s/union/viruses_union.fasta"
  [[ -s "$V" ]] && awk -v s="$s" '/^>/{print ">virus|" s "|" substr($1,2); next} {print}' "$V" >> "$TG"
done

printf 'mag\tspacer\ttarget\ttarget_type\tmismatches\tidentity\n' > "$TSV"
if [[ "$NSP" -eq 0 || ! -s "$TG" ]]; then
  echo "[crispr] $NSP spacers, $(grep -c '^>' "$TG" || true) targets: no search"; exit 0
fi
makeblastdb -in "$TG" -dbtype nucl -out "$OUT/targets" >/dev/null
blastn -task blastn-short -query "$SP" -db "$OUT/targets" -num_threads "$CPUS" -evalue 1 \
  -outfmt '6 qseqid sseqid pident length mismatch gaps qlen' -out "$OUT/hits.m6"
# full-length spacer hits with <= MAXMM mismatches + gaps
awk -F'\t' -v mm="$MAXMM" 'BEGIN{OFS="\t"} $4==$7 && ($5+$6)<=mm {
    split($1, q, "|"); t=$2; type=t; sub(/\|.*/, "", type); sub(/^[^|]*\|/, "", t)
    print q[1], $1, t, type, $5+$6, $3 }' "$OUT/hits.m6" | sort -u >> "$TSV"
echo "[crispr] $NSP spacers from MAGs, $(( $(wc -l < "$TSV") - 1 )) spacer-target hits -> $TSV"
