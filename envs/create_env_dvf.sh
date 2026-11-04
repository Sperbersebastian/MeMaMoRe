#!/usr/bin/env bash
# Create dvf_env for DeepVirFinder (not on bioconda): the old Keras 2.2.4 /
# Theano stack it was trained with, plus a shallow clone of the tool (models
# included, ~0.1 GB) into external/DeepVirFinder.
set -euo pipefail
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
micromamba create -y -n dvf_env -c conda-forge -c bioconda \
  python=3.6 numpy theano=1.0.3 keras=2.2.4 scikit-learn biopython h5py=2.10.0
DVF="$ROOT_DIR/external/DeepVirFinder"
[[ -d "$DVF/.git" ]] || git clone --depth 1 https://github.com/jessieren/DeepVirFinder.git "$DVF"
echo "[env] dvf_env created; DeepVirFinder in $DVF"
