#!/usr/bin/env bash
# Create mge_env: ISEScan (IS elements) and IntegronFinder 2 (integrons, attC/attI).
# Both ship their HMM models; no separate DB download.
set -euo pipefail
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
micromamba create -y -n mge_env -c conda-forge -c bioconda "isescan>=1.7" "integron_finder>=2.0"
echo "[env] mge_env created"
