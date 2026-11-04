#!/usr/bin/env bash
# modules/args/abricate/exec.sh
# ABRicate screens all contigs once with its bundled databases:
#   vfdb          virulence factors      -> summary/virulence_summary.tsv
#   resfinder     acquired ARGs          -> summary/args_summary.tsv (source abricate:resfinder)
#   plasmidfinder plasmid replicons      -> summary/plasmid_replicons.tsv
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
# Env:   ABRICATE_DBS (default "vfdb resfinder plasmidfinder"),
#        ABRICATE_MINID (default 80), ABRICATE_MINCOV (default 60)
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

OUT="$ROOT/SRA/args/$SAMPLE/abricate"; mkdir -p "$OUT"
DONE="$OUT/.done"
[[ "$FORCE" == "1" ]] && rm -f "$DONE"
[[ -s "$DONE" ]] && { echo "[abricate] cached"; exit 0; }

command -v abricate >/dev/null 2>&1 || { echo "[abricate] abricate not in PATH"; exit 3; }

DBS="${ABRICATE_DBS:-vfdb resfinder plasmidfinder}"
MINID="${ABRICATE_MINID:-80}"
MINCOV="${ABRICATE_MINCOV:-60}"

# ---- helper: one FASTA against every DB ----
_run_abricate() {
  local fasta="$1" context="$2" label="$3" db
  [[ -s "$fasta" ]] || { echo "[abricate] skip $context/$label (empty/missing)"; return 0; }
  for db in $DBS; do
    local out="$OUT/$db/$context/${label}.tsv"
    mkdir -p "$(dirname "$out")"
    [[ -s "$out" && "$FORCE" != "1" ]] && continue
    abricate --db "$db" --minid "$MINID" --mincov "$MINCOV" --threads "$CPUS" --nopath --quiet \
      "$fasta" > "$out.tmp" && mv -f "$out.tmp" "$out"
    echo "[abricate] $db $context/$label: $(tail -n +2 "$out" | wc -l) hits"
  done
}

total=0

# ---- all contigs (context assigned later from the contig ID) ----
CONTIGS="$ROOT/SRA/assemblies/spades/$SAMPLE/contigs.fasta"
if [[ -s "$CONTIGS" ]]; then
  _run_abricate "$CONTIGS" "contigs" "contigs"; ((total++)) || true
fi

# MAG / plasmid / virus context is assigned from the contig ID (scripts/arg_context.py)

echo "[abricate] scanned $total input(s) with: $DBS"
date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
echo "[abricate] done -> $OUT"
