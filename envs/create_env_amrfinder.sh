#!/usr/bin/env bash
# Create amrfinder_env: NCBI AMRFinderPlus (+ its DB, ~0.2 GB, via amrfinder -u)
# and argNorm (normalises ARG annotations of several tools to CARD ARO terms).
set -euo pipefail
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
micromamba create -y -n amrfinder_env -c conda-forge -c bioconda "ncbi-amrfinderplus>=4" argnorm
micromamba run -n amrfinder_env amrfinder -u
echo "[env] amrfinder_env created"
