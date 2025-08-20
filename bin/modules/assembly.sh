#!/usr/bin/env bash
set -euo pipefail

module_default_params(){ cat <<'YAML'
modules:
  assembly:
    metaspades:
      only_assembler: true
YAML
}

run_assembly(){
  local SAMPLES_TSV="$ROOT/config/samples.tsv"
  [[ -s "$SAMPLES_TSV" ]] || { echo "Missing $SAMPLES_TSV. Run ingest first."; exit 2; }

  tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r SID COUNTRY MGMT R1 R2 FASTA; do
    [[ -n "${SAMPLE:-}" && "$SID" != "$SAMPLE" ]] && continue
    ROOT="$ROOT" PARAMS_YAML="$PARAMS_YAML" SAMPLE="$SID" MODE="$MODE" R1="$R1" R2="$R2" FORCE="${FORCE:-0}" \
      bash "$ROOT/modules/assembly/spades/exec.sh"
    # after spades step per sample:
    ROOT="$ROOT" PARAMS_YAML="$PARAMS_YAML" SAMPLE="$SID" MODE="$MODE" R1="$R1" R2="$R2" FORCE="${FORCE:-0}" \
      bash "$ROOT/modules/assembly/contig_qc/exec.sh"
  done
}
