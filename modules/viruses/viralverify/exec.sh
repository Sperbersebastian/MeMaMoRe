#!/usr/bin/env bash
# modules/viruses/viralverify/exec.sh
# Run viralVerify on metaviralSPAdes contigs – extract VIRUS predictions
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

IN="$ROOT/SRA/assemblies/metaviralspades/$SAMPLE/contigs.fasta"
OUT="$ROOT/SRA/viruses/$SAMPLE/viralverify"; mkdir -p "$OUT"
if [[ ! -s "$IN" ]]; then
  echo "[vv-virus] no metaviralSPAdes contigs for $SAMPLE, skipping viralVerify"
  exit 0
fi

VV_HMM="${VV_HMM:-$ROOT/refdata/viralverify/pfam/Pfam-A.hmm}"
THR="${VV_THR:-7}"
EXTRA="${VV_EXTRA:-}"
[[ -f "$VV_HMM" ]] || { echo "[vv-virus] HMM db missing: $VV_HMM"; exit 2; }

DONE="$OUT/.done.hmm$(basename "$VV_HMM").thr${THR}"
[[ "$FORCE" == "1" ]] && rm -f "$DONE"
[[ -s "$DONE" ]] && { echo "[vv-virus] cached"; exit 0; }

# find viralverify binary
CMD=""
if command -v viralverify >/dev/null 2>&1; then
  CMD="viralverify"
elif [[ -x "$ROOT/external/viralVerify/bin/viralverify" ]]; then
  CMD="$ROOT/external/viralVerify/bin/viralverify"
else
  echo "[vv-virus] viralverify not found"; exit 4
fi

cmd=( "$CMD" -f "$IN" -o "$OUT" --hmm "$VV_HMM" -t "$CPUS" --thr "$THR" -p )
# shellcheck disable=SC2206  # extra args are intentionally word-split
[[ -n "$EXTRA" ]] && cmd+=( $EXTRA )

echo "[vv-virus] ${cmd[*]}"
"${cmd[@]}"

# Normalise virus outputs
# viralVerify writes Prediction_results_fasta/<input>_virus.fasta (older names kept as fallback)
for f in "$OUT"/Prediction_results_fasta/*_virus.fasta "$OUT"/predicted_viruses.fasta "$OUT"/virus_contigs.fasta "$OUT"/viruses.fasta; do
  [[ -s "$f" ]] && { ln -sf "${f#"$OUT"/}" "$OUT/viralverify_viruses.fasta"; break; }
done

date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
echo "[vv-virus] done -> $OUT"
