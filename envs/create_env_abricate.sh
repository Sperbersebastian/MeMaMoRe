#!/usr/bin/env bash
# Create abricate_env: ABRicate with its bundled databases (vfdb, resfinder, plasmidfinder, ...)
set -euo pipefail
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
micromamba create -y -n abricate_env -c conda-forge -c bioconda abricate
micromamba run -n abricate_env abricate --list
echo "[env] abricate_env created"
