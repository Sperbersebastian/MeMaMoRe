#!/usr/bin/env bash
# Create plasmer_env (Plasmer plasmid predictor; author channel iskoldt).
# DB (plasmerMainDB + customizedKraken2DB, Zenodo 7030675; ~14 GB packed,
# ~68 GB unpacked) is fetched on first use into refdata/plasmer. Plasmer needs
# ~32 GB of RAM + swap.
set -euo pipefail
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
micromamba create -y -n plasmer_env -c iskoldt -c bioconda -c conda-forge plasmer
echo "[env] plasmer_env created"
