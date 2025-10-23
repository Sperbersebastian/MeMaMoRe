#!/usr/bin/env bash
set -euo pipefail
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
eval "$(micromamba shell hook --shell=bash)"
micromamba create -y -n sim_env --channel-priority flexible -c bioconda -c conda-forge \
  python=3.10 art ncbi-datasets-cli entrez-direct unzip yq





