#!/usr/bin/env bash
set -euo pipefail
export MAMBA_ROOT_PREFIX="$HOME/micromamba"; export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
ENV=env_sra_tools
micromamba create -y -n "$ENV" -c conda-forge python=3.10 pyyaml ruamel.yaml click