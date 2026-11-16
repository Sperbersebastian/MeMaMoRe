#!/usr/bin/env bash
# Create gtdbtk_env (GTDB-Tk >= 2.8; needs reference data release232 or newer:
#   scripts/setup_gtdbtk_data.sh latest   -> refdata/gtdbtk/release<N>/, ~57 GB packed).
# Own env: GTDB-Tk pins (pydantic, numpy) otherwise hold back Snakemake in env_binning.
set -euo pipefail
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
micromamba create -y -n gtdbtk_env -c conda-forge -c bioconda "gtdbtk>=2.8"
echo "[env] gtdbtk_env created"
