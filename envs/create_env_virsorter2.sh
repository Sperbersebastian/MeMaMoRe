#!/usr/bin/env bash
# Create virsorter2_env micromamba environment
set -euo pipefail
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
micromamba create -y -n virsorter2_env -c conda-forge -c bioconda virsorter=2
echo "[env] virsorter2_env created"
