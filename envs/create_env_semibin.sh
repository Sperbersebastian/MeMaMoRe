#!/usr/bin/env bash
# Create semibin_env (SemiBin2, self-supervised contig binning; CPU PyTorch).
set -euo pipefail
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
micromamba create -y -n semibin_env -c conda-forge -c bioconda "semibin>=2.1" samtools
echo "[env] semibin_env created"
