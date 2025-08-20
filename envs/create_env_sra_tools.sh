#!/usr/bin/env bash
set -euo pipefail
export MAMBA_ROOT_PREFIX="$HOME/micromamba"; export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
micromamba create -y -n env_sra_tools -c conda-forge -c bioconda python=3.10 pyyaml
