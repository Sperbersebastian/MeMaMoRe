#!/usr/bin/env bash
# Create hostile_env (host read removal, used when modules.qc.host_removal = human).
# The human index (human-t2t-hla, Bowtie2) is fetched once into refdata/hostile.
set -euo pipefail
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
micromamba create -y -n hostile_env -c conda-forge -c bioconda "hostile>=2"
echo "[env] hostile_env created"
