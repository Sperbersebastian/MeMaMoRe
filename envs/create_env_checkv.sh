#!/usr/bin/env bash
# Create checkv_env micromamba environment
set -euo pipefail
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
micromamba create -y -n checkv_env -c conda-forge -c bioconda checkv diamond prodigal hmmer
echo "[env] checkv_env created"
