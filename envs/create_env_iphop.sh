#!/usr/bin/env bash
# Create iphop_env (several GB; the database is separate, ~340 GB:
#   micromamba run -n iphop_env iphop download --db_dir "$ROOT/refdata/iphop" --split)
set -euo pipefail
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
micromamba create -y -n iphop_env -c conda-forge -c bioconda "iphop>=1.3"
echo "[env] iphop_env created"
