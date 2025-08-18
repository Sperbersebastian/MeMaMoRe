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
  echo "[ingest] params at $PARAMS_YAML; MODE=$MODE; SAMPLE=${SAMPLE:-all}"
  # next step: wire programs modules/ingest/{sra,fastq,fasta}/exec.sh
}
