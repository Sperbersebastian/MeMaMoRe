#!/usr/bin/env bash
set -euo pipefail

# env: ROOT PARAMS_YAML SAMPLE COUNTRY MGMT R1_PATH R2_PATH (R2 may be empty)
WORKDIR="$ROOT/SRA"
RAW="$WORKDIR/raw_fastq/$SAMPLE"
TEST="$WORKDIR/testsets/$SAMPLE"
STATE="$WORKDIR/.state/ingest/fastq"
LOGS="$ROOT/logs"
EVENTS="$ROOT/api/status.jsonl"

mkdir -p "$RAW" "$TEST" "$STATE" "$LOGS"
ts(){ date -u +"%Y-%m-%dT%H:%M:%SZ"; }

# event: start
echo "{\"ts\":\"$(ts)\",\"module\":\"ingest\",\"program\":\"fastq\",\"sample\":\"$SAMPLE\",\"phase\":\"start\"}" >> "$EVENTS"

# copy (hardlink if possible)
copy(){ local src="$1" dst="$2"; [[ -z "$src" ]] && return 0; [[ -s "$dst" ]] && return 0; ln "$src" "$dst" 2>/dev/null || cp -f "$src" "$dst"; }

# normalize filenames
OUT1="$RAW/${SAMPLE}_1.fastq.gz"; OUT2="$RAW/${SAMPLE}_2.fastq.gz"
if [[ -n "${R2_PATH:-}" ]]; then
  copy "$R1_PATH" "$OUT1"
  copy "$R2_PATH" "$OUT2"
else
  OUT2=""
  OUT1="$RAW/${SAMPLE}.fastq.gz"
  copy "$R1_PATH" "$OUT1"
fi

# update samples.tsv (dedupe sample row, then append)
SAMPLES_TSV="$ROOT/config/samples.tsv"
[[ -s "$SAMPLES_TSV" ]] || echo -e "sample_id\tcountry\tmanagement\tR1\tR2\tFASTA" > "$SAMPLES_TSV"
awk -F'\t' -v s="$SAMPLE" 'NR==1 || $1!=s' "$SAMPLES_TSV" > "${SAMPLES_TSV}.tmp" && mv "${SAMPLES_TSV}.tmp" "$SAMPLES_TSV"
if [[ -n "$OUT2" ]]; then
  echo -e "$SAMPLE\t${COUNTRY:-NA}\t${MGMT:-NA}\t$OUT1\t$OUT2\t" >> "$SAMPLES_TSV"
else
  echo -e "$SAMPLE\t${COUNTRY:-NA}\t${MGMT:-NA}\t$OUT1\t\t" >> "$SAMPLES_TSV"
fi

# 10% testsets if enabled (case-insensitive)
mk=$(micromamba run -n env_sra_tools python "$ROOT/scripts/merge_params.py" --defaults "$PARAMS_YAML" --get modules.ingest.make_testsets)
mk=$(echo "$mk" | tr '[:upper:]' '[:lower:]')
if [[ "$mk" == "true" || "$mk" == "1" || "$mk" == "yes" ]]; then
  frac=$(micromamba run -n env_sra_tools python "$ROOT/scripts/merge_params.py" --defaults "$PARAMS_YAML" --get modules.ingest.test_fraction)
  seed=$(micromamba run -n env_sra_tools python "$ROOT/scripts/merge_params.py" --defaults "$PARAMS_YAML" --get modules.ingest.test_seed)
  if [[ -n "$OUT2" ]]; then
    T1="$TEST/${SAMPLE}_1.10p.fastq.gz"; T2="$TEST/${SAMPLE}_2.10p.fastq.gz"
    if [[ ! -s "$T1" || ! -s "$T2" ]]; then
      micromamba run -n env_sra_tools python "$ROOT/scripts/sample_fastq.py" \
        --in1 "$OUT1" --in2 "$OUT2" --out1 "$T1" --out2 "$T2" -p "$frac" -s "$seed" \
        >"$LOGS/testset_${SAMPLE}.log" 2>&1
    fi
  else
    T1="$TEST/${SAMPLE}.10p.fastq.gz"
    if [[ ! -s "$T1" ]]; then
      micromamba run -n env_sra_tools python "$ROOT/scripts/sample_fastq.py" \
        --in1 "$OUT1" --out1 "$T1" -p "$frac" -s "$seed" \
        >"$LOGS/testset_${SAMPLE}.log" 2>&1
    fi
  fi
fi

# event: done
if [[ -n "$OUT2" ]]; then
  echo "{\"ts\":\"$(ts)\",\"module\":\"ingest\",\"program\":\"fastq\",\"sample\":\"$SAMPLE\",\"phase\":\"done\",\"out\":{\"r1\":\"$OUT1\",\"r2\":\"$OUT2\"}}" >> "$EVENTS"
else
  echo "{\"ts\":\"$(ts)\",\"module\":\"ingest\",\"program\":\"fastq\",\"sample\":\"$SAMPLE\",\"phase\":\"done\",\"out\":{\"r1\":\"$OUT1\"}}" >> "$EVENTS"
fi
