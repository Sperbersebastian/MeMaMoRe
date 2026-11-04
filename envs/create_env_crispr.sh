#!/usr/bin/env bash
# Create crispr_env: MinCED (CRISPR arrays in MAGs) + BLAST (spacer matches).
set -euo pipefail
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
micromamba create -y -n crispr_env -c conda-forge -c bioconda minced blast
echo "[env] crispr_env created"
