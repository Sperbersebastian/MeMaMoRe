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
  echo "[qc] params at $PARAMS_YAML; MODE=$MODE; SAMPLE=${SAMPLE:-all}"
  # next step: call modules/qc/*/exec.sh with micromamba run
}
