#!/usr/bin/env bash
# Run metaPlasmidSPAdes on trimmed reads
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
set -euo pipefail
ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"


# --- trace who invoked me ---
printf '[mps] invoked by PPID=%s CMD="%s"\n' "$PPID" "$(ps -o command= -p "$PPID" | tail -n1)"
command -v pstree >/dev/null 2>&1 && { echo "[mps] pstree:"; pstree -salp $$ | sed -n '1,8p'; }

# Only run MPS if explicitly requested:
#  - ONLY=metaplasmidspades
#  - or MODULE_CONTEXT=plasmids (set by plasmids orchestrator)
if [[ "${ONLY:-}" != "metaplasmidspades" && "${MODULE_CONTEXT:-}" != "plasmids" ]]; then
  echo "[mps] guarded: skipping (set ONLY=metaplasmidspades or MODULE_CONTEXT=plasmids)"
  exit 0
fi


# collect reads from both layouts, ohne eval/CRLF-Probleme
# --- collect reads (only from qc/fastp) ---
shopt -s nullglob
R1=( "$ROOT/SRA/qc/fastp/$SAMPLE/"*trimmed_1*.fastq.gz )
R2=( "$ROOT/SRA/qc/fastp/$SAMPLE/"*trimmed_2*.fastq.gz )

shopt -u nullglob

[[ ${#R1[@]} -gt 0 ]] || { echo "[mps] no R1 reads found in qc/fastp"; exit 3; }


OUT="$ROOT/SRA/assemblies/metaplasmidspades/$SAMPLE"
CTG="$OUT/contigs.fasta"

# binary tolerant finden
BIN="metaplasmidspades.py"
command -v "$BIN" >/dev/null 2>&1 || BIN="metaPlasmidSPAdes.py"
command -v "$BIN" >/dev/null 2>&1 || { echo "[mps] metaPlasmidSPAdes not in PATH"; exit 2; }

[[ ${#R1[@]} -gt 0 ]] || { echo "[mps] no R1 reads found"; exit 3; }

# Rerun logic
if [[ "$FORCE" == "1" ]]; then rm -rf "$OUT"; fi
if [[ -s "$CTG" ]]; then echo "[mps] cached -> $CTG"; exit 0; fi
mkdir -p "$OUT"

# Build command (supports multi-lane; pairs by index)
cmd=( "$BIN" -o "$OUT" -t "$CPUS" )
if [[ ${#R2[@]} -gt 0 ]]; then
  [[ ${#R1[@]} -eq ${#R2[@]} ]] || { echo "[mps] R1/R2 count mismatch"; exit 4; }
  for i in "${!R1[@]}"; do
    cmd+=( -1 "${R1[$i]}" -2 "${R2[$i]}" )
  done
else
  for f in "${R1[@]}"; do cmd+=( -s "$f" ); done
fi

echo "[mps] ${cmd[*]}"
"${cmd[@]}"

[[ -s "$CTG" ]] || { echo "[mps] assembly failed, contigs not found: $CTG"; exit 5; }
echo "[mps] done -> $CTG"
