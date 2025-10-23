#!/usr/bin/env bash
# modules/plasmids/viralverify/exec.sh
# Usage: exec.sh ROOT SAMPLE [CPUS] [FORCE]
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

IN="$ROOT/SRA/assemblies/metaplasmidspades/$SAMPLE/contigs.fasta"
OUT="$ROOT/SRA/plasmids/$SAMPLE/viralverify"; mkdir -p "$OUT"
[[ -s "$IN" ]] || { echo "[vv] no input: $IN"; exit 1; }

VV_HMM="${VV_HMM:-$ROOT/refdata/viralverify/pfam/Pfam-A.hmm}"
VV_DB="${VV_DB:-}"
THR="${VV_THR:-7}"
EXTRA="${VV_EXTRA:-}"
[[ -f "$VV_HMM" ]] || { echo "[vv] HMM db missing: $VV_HMM"; exit 2; }

DONE="$OUT/.done.hmm$(basename "$VV_HMM").thr${THR}"
[[ "$FORCE" == "1" ]] && rm -f "$DONE"
[[ -s "$DONE" ]] && { echo "[vv] cached"; exit 0; }

# bevorzugt conda-CLI, sonst Repo-Binary
CMD=""
if command -v viralverify >/dev/null 2>&1; then
  CMD="viralverify"
elif [[ -x "$ROOT/external/viralVerify/bin/viralverify" ]]; then
  CMD="$ROOT/external/viralVerify/bin/viralverify"
else
  echo "[vv] viralverify not found (neither conda CLI nor external/bin)"; exit 4;
fi

# bauen
cmd=( "$CMD" -f "$IN" -o "$OUT" --hmm "$VV_HMM" -t "$CPUS" --thr "$THR" -p )
[[ -n "$VV_DB" ]] && cmd+=( --db "$VV_DB" )
[[ -n "$EXTRA" ]] && cmd+=( $EXTRA )

echo "[vv] ${cmd[*]}"
"${cmd[@]}"

# Outputs normalisieren
pl_out=""
for f in "$OUT"/predicted_plasmids.fasta "$OUT"/plasmid_contigs.fasta "$OUT"/plasmids.fasta; do
  [[ -s "$f" ]] && { pl_out="$f"; break; }
done
[[ -n "$pl_out" ]] && ln -sf "$(basename "$pl_out")" "$OUT/viralverify_plasmids.fasta" || \
  echo "[vv] WARNING: no plasmid fasta found in $OUT"

date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE"
echo "[vv] done -> $OUT"
