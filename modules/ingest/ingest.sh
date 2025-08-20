#!/usr/bin/env bash
set -euo pipefail

module_default_params(){ cat <<'YAML'
modules:
  ingest:
    make_testsets: true
    test_fraction: 0.10
    test_seed: 42
YAML
}

run_ingest(){
  # flags: --fastq <manifest> (TSV: sample_id country management R1 R2)
  local FASTQ_MANIFEST=""
  # parse trailing flags after main handled common ones
  for a in "$@"; do
    case "$a" in
      --fastq) shift; FASTQ_MANIFEST="${1:-}";;
      --fastq=*) FASTQ_MANIFEST="${a#--fastq=}";;
    esac
  done

  [[ -n "$FASTQ_MANIFEST" && -s "$FASTQ_MANIFEST" ]] || { echo "Provide --fastq <manifest.tsv>"; exit 2; }

  # validate header
  local hdr; hdr="$(head -n1 "$FASTQ_MANIFEST" | tr -d '\r')"
  for col in sample_id country management R1 R2; do
    grep -qw "$col" <<<"$hdr" || { echo "FASTQ manifest missing column: $col"; exit 2; }
  done

  # loop rows
  tail -n +2 "$FASTQ_MANIFEST" | tr -d '\r' | while IFS=$'\t' read -r SAMPLE COUNTRY MGMT R1 R2; do
    [[ -z "$SAMPLE" ]] && continue
    echo "[ingest:fastq] $SAMPLE"
    ROOT="$ROOT" PARAMS_YAML="$PARAMS_YAML" SAMPLE="$SAMPLE" COUNTRY="$COUNTRY" MGMT="$MGMT" R1_PATH="$R1" R2_PATH="$R2" \
      bash "$ROOT/modules/ingest/fastq/exec.sh" || {
        echo "error on $SAMPLE" >&2
        echo "{\"ts\":\"$(date -u +%FT%TZ)\",\"module\":\"ingest\",\"program\":\"fastq\",\"sample\":\"$SAMPLE\",\"phase\":\"error\"}" >> "$ROOT/api/status.jsonl"
        exit 3
      }
  done
}
