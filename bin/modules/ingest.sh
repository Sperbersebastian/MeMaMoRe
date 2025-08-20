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
  local FASTQ_MANIFEST=""
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --fastq) FASTQ_MANIFEST="$2"; shift 2;;
      --fastq=*) FASTQ_MANIFEST="${1#--fastq=}"; shift;;
      --set|--test|--resume|--force|--from|--only|--sample) shift 2||true;;
      *) shift;;
    esac
  done
  [[ -n "$FASTQ_MANIFEST" && -s "$FASTQ_MANIFEST" ]] || { echo "Provide --fastq <manifest.tsv>"; exit 2; }
  local hdr; hdr="$(head -n1 "$FASTQ_MANIFEST" | tr -d '\r')"
  for c in sample_id country management R1 R2; do grep -qw "$c" <<<"$hdr" || { echo "Missing $c"; exit 2; }; done
  tail -n +2 "$FASTQ_MANIFEST" | tr -d '\r' | while IFS=$'\t' read -r SAMPLE COUNTRY MGMT R1 R2; do
    [[ -z "$SAMPLE" ]] && continue
    ROOT="$ROOT" PARAMS_YAML="$PARAMS_YAML" SAMPLE="$SAMPLE" COUNTRY="$COUNTRY" MGMT="$MGMT" R1_PATH="$R1" R2_PATH="$R2" \
      bash "$ROOT/modules/ingest/fastq/exec.sh"
  done
}
