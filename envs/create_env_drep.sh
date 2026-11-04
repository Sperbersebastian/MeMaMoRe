#!/usr/bin/env bash
# Create drep_env: dRep v4 (skani for ANI) + CoverM for the project MAG catalog.
set -euo pipefail
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
micromamba create -y -n drep_env -c conda-forge -c bioconda "drep>=4" skani coverm
echo "[env] drep_env created"
