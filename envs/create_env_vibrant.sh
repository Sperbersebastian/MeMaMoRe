#!/usr/bin/env bash
# Create vibrant_env micromamba environment
set -euo pipefail
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
micromamba create -y -n vibrant_env -c conda-forge -c bioconda vibrant=1.2.1 prodigal
echo "[env] vibrant_env created"
