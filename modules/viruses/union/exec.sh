#!/usr/bin/env bash
# modules/viruses/union/exec.sh
# Collect + deduplicate virus FASTAs from viralVerify, geNomad, VirSorter2, VIBRANT
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

OUT="$ROOT/SRA/viruses/$SAMPLE/union"; mkdir -p "$OUT"
DONE="$OUT/.done"
[[ "$FORCE" == "1" ]] && rm -f "$DONE"
[[ -s "$DONE" ]] && { echo "[virus-union] cached"; exit 0; }

CONCAT="$OUT/viruses_concat.fasta"
: > "$CONCAT"

# Collect virus FASTAs from all tools
count=0
for src in \
  "$ROOT/SRA/viruses/$SAMPLE/viralverify/viralverify_viruses.fasta" \
  "$ROOT/SRA/viruses/$SAMPLE/genomad/viruses.fna" \
  "$ROOT/SRA/viruses/$SAMPLE/virsorter2/viruses.fasta" \
  "$ROOT/SRA/viruses/$SAMPLE/vibrant/viruses.fasta"
do
  if [[ -s "$src" ]]; then
    echo "[virus-union] adding $(basename "$(dirname "$src")"): $(grep -c '^>' "$src" || echo 0) seqs"
    cat "$src" >> "$CONCAT"
    ((count++)) || true
  else
    echo "[virus-union] skip $(basename "$(dirname "$src")") (not found)"
  fi
done

[[ $count -gt 0 ]] || { echo "[virus-union] no virus FASTAs found from any tool"; exit 1; }

TOTAL=$(grep -c '^>' "$CONCAT" || echo 0)
echo "[virus-union] total concatenated: $TOTAL seqs from $count tools"

# Deduplicate with CD-HIT at 95% identity (standard for viral genomes)
UNION="$OUT/viruses_union.fasta"
if command -v cd-hit-est >/dev/null 2>&1; then
  cd-hit-est -i "$CONCAT" -o "$UNION" -c 0.95 -aS 0.85 -M 0 -T "$CPUS" -d 0
  NUNIQ=$(grep -c '^>' "$UNION" || echo 0)
  echo "[virus-union] deduplicated: $NUNIQ / $TOTAL"
else
  echo "[virus-union] cd-hit-est not found, using concat as union"
  cp "$CONCAT" "$UNION"
fi

date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
echo "[virus-union] done -> $UNION"
