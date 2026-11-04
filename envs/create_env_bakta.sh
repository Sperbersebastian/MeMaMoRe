#!/usr/bin/env bash
# Create bakta_env (optional MAG annotation, modules.binning.bakta: true). The DB
# (light ~1.4 GB packed / 3.4 GB, full ~37 / 71 GB) is fetched on first use into
# refdata/bakta by modules/binning/bakta/exec.sh.
set -euo pipefail
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
micromamba create -y -n bakta_env -c conda-forge -c bioconda "bakta>=1.9"
echo "[env] bakta_env created"
