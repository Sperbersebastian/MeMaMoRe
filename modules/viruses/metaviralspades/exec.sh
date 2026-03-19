#!/usr/bin/env bash
# modules/viruses/metaviralspades/exec.sh
# Run metaviralSPAdes on trimmed reads for virus-enriched assembly
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

# --- collect trimmed reads from qc/fastp ---
shopt -s nullglob
R1=( "$ROOT/SRA/qc/fastp/$SAMPLE/"*trimmed_1*.fastq.gz )
R2=( "$ROOT/SRA/qc/fastp/$SAMPLE/"*trimmed_2*.fastq.gz )
shopt -u nullglob

[[ ${#R1[@]} -gt 0 ]] || { echo "[mvs] no R1 reads found in qc/fastp"; exit 3; }

OUT="$ROOT/SRA/assemblies/metaviralspades/$SAMPLE"
CTG="$OUT/contigs.fasta"

# binary lookup (SPAdes ships both names)
BIN="metaviralspades.py"
command -v "$BIN" >/dev/null 2>&1 || BIN="metaviralSPAdes.py"
command -v "$BIN" >/dev/null 2>&1 || { echo "[mvs] metaviralSPAdes not in PATH"; exit 2; }

# Rerun logic
if [[ "$FORCE" == "1" ]]; then rm -rf "$OUT"; fi
if [[ -s "$CTG" ]]; then echo "[mvs] cached -> $CTG"; exit 0; fi
mkdir -p "$OUT"

# Build command (supports multi-lane; pairs by index)
cmd=( "$BIN" -o "$OUT" -t "$CPUS" )
if [[ ${#R2[@]} -gt 0 ]]; then
  [[ ${#R1[@]} -eq ${#R2[@]} ]] || { echo "[mvs] R1/R2 count mismatch"; exit 4; }
  for i in "${!R1[@]}"; do
    cmd+=( -1 "${R1[$i]}" -2 "${R2[$i]}" )
  done
else
  for f in "${R1[@]}"; do cmd+=( -s "$f" ); done
fi

echo "[mvs] ${cmd[*]}"
"${cmd[@]}"

[[ -s "$CTG" ]] || {
  # Try scaffolds as fallback
  if [[ -s "$OUT/scaffolds.fasta" ]]; then
    echo "[mvs] contigs.fasta not found, using scaffolds.fasta"
    ln -sf scaffolds.fasta "$CTG"
  elif [[ -s "$OUT/before_rr.fasta" ]]; then
    echo "[mvs] contigs.fasta not found, using before_rr.fasta"
    ln -sf before_rr.fasta "$CTG"
  else
    echo "[mvs] WARNING: no viral contigs assembled for $SAMPLE (this is normal for some samples)"
    exit 0
  fi
}
echo "[mvs] done -> $CTG"
