#!/usr/bin/env bash
set -euo pipefail
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
micromamba create -y -n env_mapping_coverm -c bioconda -c conda-forge bwa-mem2 samtools coverm
