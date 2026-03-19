#!/usr/bin/env bash
# modules/viruses/virsorter2/exec.sh
# Run VirSorter2 on metaSPAdes contigs
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

# prefer metaSPAdes contigs
IN="$ROOT/SRA/assemblies/spades/$SAMPLE/contigs.fasta"
[[ -s "$IN" ]] || IN="$ROOT/SRA/assemblies/metaviralspades/$SAMPLE/contigs.fasta"
[[ -s "$IN" ]] || { echo "[vs2] no contigs found for $SAMPLE"; exit 1; }

OUT="$ROOT/SRA/viruses/$SAMPLE/virsorter2"; mkdir -p "$OUT"
DONE="$OUT/.done"
[[ "$FORCE" == "1" ]] && rm -f "$DONE"
[[ -s "$DONE" ]] && { echo "[vs2] cached"; exit 0; }

command -v virsorter >/dev/null 2>&1 || { echo "[vs2] virsorter not in PATH"; exit 3; }

# Database setup – auto-download if missing
VS2_DB="${VS2_DB:-$ROOT/refdata/virsorter2/db}"
if [[ ! -d "$VS2_DB" || ! -f "$VS2_DB/Done_all_setup" ]]; then
  echo "[vs2] DB missing -> downloading to $VS2_DB"
  mkdir -p "$VS2_DB"
  virsorter setup -d "$VS2_DB" -j "$CPUS"
  [[ -f "$VS2_DB/Done_all_setup" ]] || { echo "[vs2] DB setup failed"; exit 2; }
  echo "[vs2] DB ready: $VS2_DB"
fi

# VirSorter2 run (all groups, include-groups recommended for metagenomes)
cmd=( virsorter run
  -i "$IN"
  -w "$OUT"
  -d "$VS2_DB"
  --include-groups "dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae"
  -j "$CPUS"
  --min-length 1000
  all
)
echo "[vs2] ${cmd[*]}"
"${cmd[@]}"

# normalise outputs
VS2_FA="$OUT/final-viral-combined.fa"
if [[ -s "$VS2_FA" ]]; then
  ln -sf "final-viral-combined.fa" "$OUT/viruses.fasta"
else
  echo "[vs2] WARNING: no viral sequences found"
fi

date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
echo "[vs2] done -> $OUT"
