#!/usr/bin/env bash
# Create gunc_env (GUNC chimerism check for MAGs). The reference DB is fetched
# on first use by modules/binning/gunc/exec.sh into refdata/gunc (ProGenomes 2.1).
set -euo pipefail
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
micromamba create -y -n gunc_env -c conda-forge -c bioconda "gunc>=1.0.6"
echo "[env] gunc_env created"
