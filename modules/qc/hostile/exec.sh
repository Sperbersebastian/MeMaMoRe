#!/usr/bin/env bash
# modules/qc/hostile/exec.sh
# Host read removal with Hostile after fastp (modules.qc.host_removal = human).
# Index human-t2t-hla (T2T-CHM13v2.0 + IPD-IMGT/HLA, Bowtie2) is fetched once into
# refdata/hostile (HOSTILE_CACHE_DIR). Output: SRA/qc/hostile/<s>/<s>_clean_{1,2}.fastq.gz
# + hostile.json (reads in/out, host fraction). Downstream steps pick these up via
# bin/lib/reads.sh.
# env in: ROOT PARAMS_YAML SAMPLE IN1 [IN2] [FORCE]
set -euo pipefail

: "${ROOT:?}"; : "${SAMPLE:?}"; : "${IN1:?}"; IN2="${IN2:-}"
OUT="$ROOT/SRA/qc/hostile/$SAMPLE"
LOG="$ROOT/logs/qc_hostile_${SAMPLE}.log"
mkdir -p "$OUT" "$ROOT/logs"
# shellcheck source=bin/lib/cfg.sh
source "$ROOT/bin/lib/cfg.sh"
THREADS="$(cfg global.threads 8)"
INDEX="$(cfg modules.qc.host_index human-t2t-hla)"
export HOSTILE_CACHE_DIR="${HOSTILE_CACHE_DIR:-$ROOT/refdata/hostile}"
mkdir -p "$HOSTILE_CACHE_DIR/tmp"
# the first run downloads the index (~4 GB) via a temp file: keep it off the system disk
export TMPDIR="$HOSTILE_CACHE_DIR/tmp"

C1="$OUT/${SAMPLE}_clean_1.fastq.gz"; C2="$OUT/${SAMPLE}_clean_2.fastq.gz"
if [[ -s "$C1" && "$C1" -nt "$IN1" && "${FORCE:-0}" != "1" ]]; then
  echo "[hostile] cached -> $OUT"; exit 0
fi

micromamba env list --json | grep -q '/hostile_env"' || bash "$ROOT/envs/create_env_hostile.sh"

TMP="$OUT/tmp"; rm -rf "$TMP"; mkdir -p "$TMP"
cmd=(hostile clean --fastq1 "$IN1" --index "$INDEX" --aligner bowtie2 -t "$THREADS" -o "$TMP" --force)
[[ -n "$IN2" ]] && cmd+=(--fastq2 "$IN2")
echo "[hostile] ${cmd[*]}" >>"$LOG"
micromamba run -n hostile_env "${cmd[@]}" > "$OUT/hostile.json" 2>>"$LOG"

# fixed names for downstream steps (Hostile appends .clean_1/.clean_2 to the input names)
shopt -s nullglob
o1=( "$TMP"/*.clean_1.fastq.gz ); o2=( "$TMP"/*.clean_2.fastq.gz )
shopt -u nullglob
[[ ${#o1[@]} -eq 1 ]] || { echo "[hostile] ERROR: no cleaned R1 in $TMP (see $LOG)"; exit 3; }
mv -f "${o1[0]}" "$C1"
if [[ -n "$IN2" ]]; then
  [[ ${#o2[@]} -eq 1 ]] || { echo "[hostile] ERROR: no cleaned R2 in $TMP"; exit 3; }
  mv -f "${o2[0]}" "$C2"
fi
rm -rf "$TMP"
python3 - "$OUT/hostile.json" <<'PY' || true
import json, sys
d = json.load(open(sys.argv[1]))
d = d[0] if isinstance(d, list) else d
print(f"[hostile] reads in {d.get('reads_in')}, out {d.get('reads_out')}, "
      f"host removed {d.get('reads_removed')} ({d.get('reads_removed_proportion')})")
PY
