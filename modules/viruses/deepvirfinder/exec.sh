#!/usr/bin/env bash
# modules/viruses/deepvirfinder/exec.sh
# DeepVirFinder (5th virus predictor) on the metaSPAdes contigs >= modules.viruses.min_len.
# Contigs with score >= modules.viruses.dvf_min_score (0.9) and p < dvf_max_p
# (0.05) go to deepvirfinder/viruses.fasta, which the virus union picks up.
# Tool: external/DeepVirFinder (git clone by envs/create_env_dvf.sh), dvf_env.
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"
OUT="$ROOT/SRA/viruses/$SAMPLE/deepvirfinder"; mkdir -p "$OUT"
DONE="$OUT/.done"
[[ "$FORCE" == "1" ]] && rm -f "$DONE"
[[ -s "$DONE" ]] && { echo "[dvf] cached"; exit 0; }

# shellcheck source=bin/lib/cfg.sh
source "$ROOT/bin/lib/cfg.sh"
MINLEN="$(cfg modules.viruses.min_len 1000)"
SCORE="$(cfg modules.viruses.dvf_min_score 0.9)"
MAXP="$(cfg modules.viruses.dvf_max_p 0.05)"
DVF="${DVF_DIR:-$ROOT/external/DeepVirFinder}"
[[ -s "$DVF/dvf.py" ]] || { echo "[dvf] $DVF/dvf.py missing (run envs/create_env_dvf.sh)"; exit 3; }

IN="$ROOT/SRA/assemblies/spades/$SAMPLE/contigs.fasta"
[[ -s "$IN" ]] || { echo "[dvf] no contigs for $SAMPLE"; exit 1; }

# Theano compiles into ~/.theano by default; keep that next to the results
export THEANO_FLAGS="${THEANO_FLAGS:-base_compiledir=$OUT/.theano}"
python "$DVF/dvf.py" -i "$IN" -o "$OUT/raw" -l "$MINLEN" -c "$CPUS" > "$OUT/dvf.log" 2>&1
PRED="$(find "$OUT/raw" -name '*dvfpred.txt' -print -quit)"
[[ -s "$PRED" ]] || { echo "[dvf] ERROR: no *dvfpred.txt (see $OUT/dvf.log)"; exit 4; }

# columns: name len score pvalue
awk -F'\t' -v s="$SCORE" -v p="$MAXP" 'NR>1 && $3>=s && $4<p {print $1}' "$PRED" > "$OUT/viral.ids"
awk 'NR==FNR{keep[$1]=1; next} /^>/{id=substr($1,2); k=(id in keep)} k' "$OUT/viral.ids" "$IN" > "$OUT/viruses.fasta"
cp -f "$PRED" "$OUT/dvfpred.txt"
echo "[dvf] $(wc -l < "$OUT/viral.ids") contigs with score >= $SCORE and p < $MAXP -> $OUT/viruses.fasta"
date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
