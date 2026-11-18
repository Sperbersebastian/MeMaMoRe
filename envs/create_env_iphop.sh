#!/usr/bin/env bash
# Create iphop_env (several GB; the database is separate, ~340 GB:
#   scripts/setup_iphop_db.sh  (chunk-wise, needs only the unpacked size on disk))
set -euo pipefail
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
micromamba create -y -n iphop_env -c conda-forge -c bioconda "iphop>=1.3"
echo "[env] iphop_env created"
