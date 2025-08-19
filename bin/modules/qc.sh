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
  local SAMPLES_TSV="$ROOT/config/samples.tsv"
  [[ -s "$SAMPLES_TSV" ]] || { echo "Missing $SAMPLES_TSV. Run ingest first."; exit 2; }

  # loop samples
  tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r SID COUNTRY MGMT R1 R2 FASTA; do
    [[ -n "${SAMPLE:-}" && "$SID" != "$SAMPLE" ]] && continue

    # choose inputs (prefer testsets in --test)
    if [[ "${MODE:-full}" == "test" ]]; then
      if [[ -n "$R2" ]]; then
        IN1="$ROOT/SRA/testsets/$SID/${SID}_1.10p.fastq.gz"
        IN2="$ROOT/SRA/testsets/$SID/${SID}_2.10p.fastq.gz"
      else
        IN1="$ROOT/SRA/testsets/$SID/${SID}.10p.fastq.gz"
        IN2=""
      fi
      if [[ ! -s "$IN1" ]]; then
        echo "[qc] WARNING: testset missing for $SID → falling back to raw"
        IN1="$R1"; IN2="$R2"
      fi
    else
      IN1="$R1"; IN2="$R2"
    fi
    echo "[qc] $SID MODE=$MODE IN1=$IN1 IN2=$IN2"

    # 1) fastp
    ROOT="$ROOT" PARAMS_YAML="$PARAMS_YAML" SAMPLE="$SID" MODE="$MODE" IN1="$IN1" IN2="$IN2" FORCE="${FORCE:-0}" \
      bash "$ROOT/modules/qc/fastp/exec.sh"

    # 2) fastqc on trimmed outputs
    if [[ -n "$IN2" ]]; then
      T1="$ROOT/SRA/qc/fastp/$SID/${SID}_trimmed_1.fastq.gz"
      T2="$ROOT/SRA/qc/fastp/$SID/${SID}_trimmed_2.fastq.gz"
    else
      T1="$ROOT/SRA/qc/fastp/$SID/${SID}_trimmed_1.fastq.gz"
      T2=""
    fi
    ROOT="$ROOT" PARAMS_YAML="$PARAMS_YAML" SAMPLE="$SID" IN1="$T1" IN2="$T2" \
      bash "$ROOT/modules/qc/fastqc/exec.sh"
  done

  # 3) multiqc once at the end
  ROOT="$ROOT" PARAMS_YAML="$PARAMS_YAML" bash "$ROOT/modules/qc/multiqc/exec.sh"
}
