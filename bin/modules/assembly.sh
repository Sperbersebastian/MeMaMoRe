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
  # Use the canonical 3-col manifest produced by ingest
  local MANI="$ROOT/SRA/reads/manifest.tsv"
  [[ -s "$MANI" ]] || { echo "Missing $MANI. Run ingest first."; exit 2; }

  # sample \t fq1 \t fq2
  tail -n +2 "$MANI" | tr -d '\r' | while IFS=$'\t' read -r SID R1 R2; do
    [[ -n "${SAMPLE:-}" && "$SID" != "$SAMPLE" ]] && continue

    # choose inputs (prefer testsets in --test if present)
    local IN1 IN2
    if [[ "${MODE:-full}" == "test" ]]; then
      if [[ -n "${R2:-}" ]]; then
        IN1="$ROOT/SRA/testsets/$SID/${SID}_1.10p.fastq.gz"
        IN2="$ROOT/SRA/testsets/$SID/${SID}_2.10p.fastq.gz"
      else
        IN1="$ROOT/SRA/testsets/$SID/${SID}.10p.fastq.gz"
        IN2=""
      fi
      # fallback to raw if test subset missing
      [[ -s "$IN1" ]] || { IN1="$R1"; IN2="${R2:-}"; }
    else
      IN1="$R1"; IN2="${R2:-}"
    fi

    echo "[assembly] $SID MODE=$MODE R1=$IN1 R2=${IN2:-}"

    # 1) SPAdes/metaspades per sample
    ROOT="$ROOT" PARAMS_YAML="$PARAMS_YAML" SAMPLE="$SID" MODE="$MODE" R1="$IN1" R2="${IN2:-}" FORCE="${FORCE:-0}" \
      bash "$ROOT/modules/assembly/spades/exec.sh"

    # 2) Contig QC per sample (expects outputs from step 1)
    ROOT="$ROOT" PARAMS_YAML="$PARAMS_YAML" SAMPLE="$SID" MODE="$MODE" R1="$IN1" R2="${IN2:-}" FORCE="${FORCE:-0}" \
      bash "$ROOT/modules/assembly/contig_qc/exec.sh"
  done
}
