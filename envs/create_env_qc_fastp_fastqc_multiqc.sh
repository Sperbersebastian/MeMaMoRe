#!/usr/bin/env bash
set -euo pipefail
export MAMBA_ROOT_PREFIX="$HOME/micromamba"; export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
micromamba create -y -n env_qc_fastp_fastqc_multiqc -c bioconda -c conda-forge fastp fastqc multiqc
