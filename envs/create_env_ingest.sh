#!/usr/bin/env bash
set -euo pipefail
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
ENV=env_ingest
# sra-tools + pigz are needed for SRR downloads (fasterq-dump, compression)
micromamba create -y -n "$ENV" -c conda-forge -c bioconda python=3.10 pandas pyyaml click sra-tools pigz
