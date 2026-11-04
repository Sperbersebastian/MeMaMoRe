#!/usr/bin/env bash
# modules/viruses/vibrant/exec.sh
# Run VIBRANT on metaSPAdes contigs
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

# prefer metaSPAdes contigs
IN="$ROOT/SRA/assemblies/spades/$SAMPLE/contigs.fasta"
[[ -s "$IN" ]] || IN="$ROOT/SRA/assemblies/metaviralspades/$SAMPLE/contigs.fasta"
[[ -s "$IN" ]] || { echo "[vibrant] no contigs found for $SAMPLE"; exit 1; }

OUT="$ROOT/SRA/viruses/$SAMPLE/vibrant"; mkdir -p "$OUT"
DONE="$OUT/.done"
[[ "$FORCE" == "1" ]] && rm -f "$DONE"
[[ -s "$DONE" ]] && { echo "[vibrant] cached"; exit 0; }

command -v VIBRANT_run.py >/dev/null 2>&1 || { echo "[vibrant] VIBRANT_run.py not in PATH"; exit 3; }

# VIBRANT requires its HMM database to be setup before first run.
# Override with VIBRANT_DB; defaults to the project refdata dir.
VIBRANT_DB="${VIBRANT_DB:-$ROOT/refdata/vibrant/db}"
export VIBRANT_DATA_PATH="$VIBRANT_DB"
DBD="$VIBRANT_DB/databases"
HMMS=(KEGG_profiles_prokaryotes.HMM Pfam-A_v32.HMM VOGDB94_phage.HMM)
_have_hmm(){ [[ -s "$DBD/$1" || -s "$DBD/$1.h3p" ]]; }
if ! _have_hmm "${HMMS[0]}" || ! _have_hmm "${HMMS[1]}" || ! _have_hmm "${HMMS[2]}" || [[ ! -d "$VIBRANT_DB/files" ]]; then
  echo "[vibrant] DB missing or incomplete -> downloading and setting up"
  command -v download-db.sh >/dev/null 2>&1 || { echo "[vibrant] download-db.sh not found"; exit 4; }
  download-db.sh "$VIBRANT_DB" || { echo "[vibrant] DB setup failed"; exit 5; }
fi
# an interrupted setup can leave profiles un-pressed ("could not identify KEGG HMM files")
for h in "${HMMS[@]}"; do
  if [[ -s "$DBD/$h" ]] && { [[ ! -s "$DBD/$h.h3f" ]] || [[ ! -s "$DBD/$h.h3i" ]] || [[ ! -s "$DBD/$h.h3m" ]] || [[ ! -s "$DBD/$h.h3p" ]]; }; then
    echo "[vibrant] hmmpress $h (completing DB setup)"
    hmmpress -f "$DBD/$h" >/dev/null
  fi
done

# VIBRANT expects the input basename for output naming
IN_BASE="$(basename "$IN" .fasta)"

# shellcheck source=bin/lib/cfg.sh
source "$ROOT/bin/lib/cfg.sh"
MINLEN="$(cfg modules.viruses.min_len 1000)"
cmd=( VIBRANT_run.py
  -i "$IN"
  -folder "$OUT"
  -d "$DBD"
  -m "$VIBRANT_DB/files"
  -t "$CPUS"
  -l "$MINLEN"
  -o "$(cfg modules.viruses.vibrant_min_orfs 4)"
)
case "$(cfg modules.viruses.vibrant_virome false)" in true|True|1|yes) cmd+=( -virome ) ;; esac
echo "[vibrant] ${cmd[*]}"
"${cmd[@]}"
# VIBRANT reports setup/input problems in its log but still exits 0
VLOG="$OUT/VIBRANT_log_run_${IN_BASE}.log"
if grep -q "VIBRANT error" "$VLOG" 2>/dev/null; then
  echo "[vibrant] $(grep -m1 "VIBRANT error" "$VLOG")"; exit 6
fi

# normalise outputs – VIBRANT creates a nested structure
VIBRANT_DIR="$OUT/VIBRANT_${IN_BASE}"
PHAGE_FA="$VIBRANT_DIR/VIBRANT_phages_${IN_BASE}/${IN_BASE}.phages_combined.fna"
if [[ -s "$PHAGE_FA" ]]; then
  ln -sf "$PHAGE_FA" "$OUT/viruses.fasta"
else
  # try alternative output path
  ALT_FA="$(find "$OUT" -name '*.phages_combined.fna' -print -quit 2>/dev/null || true)"
  if [[ -s "${ALT_FA:-}" ]]; then
    ln -sf "$ALT_FA" "$OUT/viruses.fasta"
  else
    echo "[vibrant] WARNING: no viral sequences found"
  fi
fi

date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
echo "[vibrant] done -> $OUT"
