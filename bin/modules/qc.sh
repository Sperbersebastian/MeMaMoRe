#!/usr/bin/env bash
set -euo pipefail

module_default_params(){ cat <<'YAML'
modules:
  qc:
    fastp:
      length_required: 50
      qualified_quality_phred: 15
      unqualified_percent_limit: 40
YAML
}

run_qc(){
  local MANI="$ROOT/SRA/reads/manifest.tsv"
  [[ -s "$MANI" ]] || { echo "Missing $MANI. Run ingest first."; exit 2; }

  tail -n +2 "$MANI" | tr -d '\r' | while IFS=$'\t' read -r SID R1 R2; do
    [[ -n "${SAMPLE:-}" && "$SID" != "$SAMPLE" ]] && continue

    # choose inputs
    if [[ "${MODE:-full}" == "test" ]]; then
      if [[ -n "$R2" ]]; then
        IN1="$ROOT/SRA/testsets/$SID/${SID}_1.10p.fastq.gz"
        IN2="$ROOT/SRA/testsets/$SID/${SID}_2.10p.fastq.gz"
      else
        IN1="$ROOT/SRA/testsets/$SID/${SID}.10p.fastq.gz"
        IN2=""
      fi
      [[ -s "$IN1" ]] || { IN1="$R1"; IN2="$R2"; }
    else
      IN1="$R1"; IN2="$R2"
    fi

    echo "[qc] $SID MODE=$MODE IN1=$IN1 IN2=${IN2:-}"

    ROOT="$ROOT" PARAMS_YAML="$PARAMS_YAML" SAMPLE="$SID" MODE="$MODE" IN1="$IN1" IN2="${IN2:-}" FORCE="${FORCE:-0}" \
      bash "$ROOT/modules/qc/fastp/exec.sh"

    # trimmed outputs
    T1="$ROOT/SRA/qc/fastp/$SID/${SID}_trimmed_1.fastq.gz"
    T2=""
    [[ -n "${IN2:-}" ]] && T2="$ROOT/SRA/qc/fastp/$SID/${SID}_trimmed_2.fastq.gz"

    ROOT="$ROOT" PARAMS_YAML="$PARAMS_YAML" SAMPLE="$SID" IN1="$T1" IN2="${T2:-}" \
      bash "$ROOT/modules/qc/fastqc/exec.sh"
  done

  ROOT="$ROOT" PARAMS_YAML="$PARAMS_YAML" bash "$ROOT/modules/qc/multiqc/exec.sh"
}

