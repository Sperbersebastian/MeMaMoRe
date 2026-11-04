#!/usr/bin/env bash
# modules/assembly/genes/exec.sh
# Central gene prediction (Pyrodigal, metagenome mode) on the assembly contigs
# >= global.contig_min_len. Shared by DeepARG (genes.ffn), RGI and MAGScoT
# (genes.faa) and later gene-based steps; geNomad, VIBRANT, VirSorter2, CheckV,
# CheckM2 and GTDB-Tk keep their own gene calling.
# Outputs SRA/assemblies/genes/<sample>/genes.{faa,ffn,gff}; gene IDs <contig>_<n>.
# env in: ROOT PARAMS_YAML SAMPLE [FORCE]
set -euo pipefail

: "${ROOT:?}"; : "${SAMPLE:?}"; : "${PARAMS_YAML:?}"
OUTDIR="$ROOT/SRA/assemblies/genes/$SAMPLE"
LOG="$ROOT/logs/assembly_genes_${SAMPLE}.log"
EVENTS="$ROOT/api/status.jsonl"
mkdir -p "$OUTDIR" "$ROOT/logs" "$(dirname "$EVENTS")"

ts(){ date -u +"%Y-%m-%dT%H:%M:%SZ"; }
getp(){ micromamba run -n env_sra_tools python "$ROOT/scripts/merge_params.py" --defaults "$PARAMS_YAML" --get "$1"; }
THREADS="$(getp global.threads)"
MINLEN="$(getp global.contig_min_len)"

ASM="$ROOT/SRA/assemblies/spades/$SAMPLE"
IN="$ASM/contigs.len${MINLEN}.fasta"
[[ -s "$IN" ]] || { echo "[genes] no contigs >= ${MINLEN} bp for $SAMPLE ($IN)" >&2; exit 2; }

FAA="$OUTDIR/genes.faa"; FFN="$OUTDIR/genes.ffn"; GFF="$OUTDIR/genes.gff"
if [[ -s "$FAA" && -s "$FFN" && -s "$GFF" && "${FORCE:-0}" != "1" && "$FAA" -nt "$IN" ]]; then
  echo "[genes] cached -> $OUTDIR"; exit 0
fi

echo "{\"ts\":\"$(ts)\",\"module\":\"assembly\",\"program\":\"pyrodigal\",\"sample\":\"$SAMPLE\",\"phase\":\"start\"}" >> "$EVENTS"
micromamba run -n env_assembly_core pyrodigal -p meta -i "$IN" \
  -a "$FAA.tmp" -d "$FFN.tmp" -f gff -o "$GFF.tmp" -j "$THREADS" >"$LOG" 2>&1
mv -f "$FAA.tmp" "$FAA"; mv -f "$FFN.tmp" "$FFN"; mv -f "$GFF.tmp" "$GFF"

N=$(grep -c '^>' "$FAA" || true)
echo "[genes] $N genes on contigs >= ${MINLEN} bp -> $OUTDIR"
echo "{\"ts\":\"$(ts)\",\"module\":\"assembly\",\"program\":\"pyrodigal\",\"sample\":\"$SAMPLE\",\"phase\":\"done\",\"out\":{\"faa\":\"$FAA\",\"genes\":$N}}" >> "$EVENTS"
