#!/usr/bin/env bash
set -euo pipefail
# env in: ROOT PARAMS_YAML SAMPLE MODE R1 R2 FORCE

OUTDIR="$ROOT/SRA/assemblies/spades/$SAMPLE"
LOGS="$ROOT/logs"
EVENTS="$ROOT/api/status.jsonl"
mkdir -p "$OUTDIR" "$LOGS"
[[ "${FORCE:-0}" == "1" ]] && rm -rf "$OUTDIR"/*

ts(){ date -u +"%Y-%m-%dT%H:%M:%SZ"; }
echo "{\"ts\":\"$(ts)\",\"module\":\"assembly\",\"program\":\"spades\",\"sample\":\"$SAMPLE\",\"phase\":\"start\"}" >> "$EVENTS"

getp(){ micromamba run -n env_sra_tools python "$ROOT/scripts/merge_params.py" --defaults "$PARAMS_YAML" --get "$1"; }
THREADS=$(getp global.threads)
MINLEN=$(getp global.contig_min_len)

LOG="$LOGS/assembly_spades_${SAMPLE}.log"
CONTIGS="$OUTDIR/contigs.fasta"
FILTERED="$OUTDIR/contigs.len${MINLEN}.fasta"

# ── Skip logic ────────────────────────────────────────────────────────────
if [[ -s "$FILTERED" && "${FORCE:-0}" != "1" ]]; then
  echo "[assembly] Reusing existing filtered contigs for $SAMPLE: $FILTERED" >>"$LOG"
  echo "{\"ts\":\"$(ts)\",\"module\":\"assembly\",\"program\":\"spades\",\"sample\":\"$SAMPLE\",\"phase\":\"done\",\"out\":{\"contigs\":\"$FILTERED\"},\"reuse\":true}" >> "$EVENTS"
  exit 0
fi

# choose inputs
IN1=""; IN2=""
if [[ "${MODE:-full}" == "test" ]]; then
  if [[ -n "${R2:-}" ]]; then
    IN1="$ROOT/SRA/testsets/$SAMPLE/${SAMPLE}_1.10p.fastq.gz"
    IN2="$ROOT/SRA/testsets/$SAMPLE/${SAMPLE}_2.10p.fastq.gz"
  else
    IN1="$ROOT/SRA/testsets/$SAMPLE/${SAMPLE}.10p.fastq.gz"
  fi
else
  T1="$ROOT/SRA/qc/fastp/$SAMPLE/${SAMPLE}_trimmed_1.fastq.gz"
  T2="$ROOT/SRA/qc/fastp/$SAMPLE/${SAMPLE}_trimmed_2.fastq.gz"
  if [[ -s "$T1" ]]; then IN1="$T1"; [[ -s "$T2" ]] && IN2="$T2"; else IN1="$R1"; IN2="${R2:-}"; fi
fi

# ── Run SPAdes only if needed ─────────────────────────────────────────────
if [[ -s "$CONTIGS" && "${FORCE:-0}" != "1" ]]; then
  echo "[assembly] Reusing existing SPAdes contigs for $SAMPLE: $CONTIGS" >>"$LOG"
else
  if [[ -n "${FORCE:-}" && "$FORCE" == "1" ]]; then
    rm -f "$OUTDIR"/*
  fi
  if [[ -n "$IN2" ]]; then
    micromamba run -n env_assembly_spades metaspades.py -1 "$IN1" -2 "$IN2" -o "$OUTDIR" -t "$THREADS" --only-assembler >"$LOG" 2>&1
  else
    micromamba run -n env_assembly_spades metaspades.py -s "$IN1" -o "$OUTDIR" -t "$THREADS" --only-assembler >"$LOG" 2>&1
  fi
fi

# ── Filter contigs by length (skip if already done unless FORCE) ──────────
if [[ -s "$CONTIGS" ]]; then
  if [[ -s "$FILTERED" && "${FORCE:-0}" != "1" ]]; then
    echo "[assembly] Filtered file already exists, skipping: $FILTERED" >>"$LOG"
  else
    export OUTDIR MINLEN
    micromamba run -n env_sra_tools python - <<'PY'
import os
outdir = os.environ["OUTDIR"]
minlen = int(os.environ.get("MINLEN","1000"))
path = os.path.join(outdir, "contigs.fasta")
outpath = os.path.join(outdir, f"contigs.len{minlen}.fasta")

def write_wrapped(fw, sid, seq, w=80):
    fw.write(f">{sid}\n")
    for i in range(0, len(seq), w):
        fw.write(seq[i:i+w] + "\n")

sid=None; buf=[]
with open(path) as fr, open(outpath, "w") as fw:
    def flush():
        global sid, buf
        if sid is None: return
        seq="".join(buf)
        if len(seq)>=minlen:
            write_wrapped(fw, sid, seq)
    for line in fr:
        if line.startswith(">"):
            flush(); sid=line[1:].strip(); buf=[]
        else:
            buf.append(line.strip())
    flush()
PY
  fi
fi

echo "{\"ts\":\"$(ts)\",\"module\":\"assembly\",\"program\":\"spades\",\"sample\":\"$SAMPLE\",\"phase\":\"done\",\"out\":{\"contigs\":\"$FILTERED\"}}" >> "$EVENTS"
