#!/usr/bin/env bash
set -euo pipefail
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
eval "$(micromamba shell hook --shell=bash)"

ENV_NAME="env_binning"
PKGS=(metabat2 snakemake yq checkm2 gtdbtk barrnap tRNAscan-SE samtools seqkit conda conda-libmamba-solver r-base r-optparse r-dplyr r-readr r-funr r-digest hmmer prodigal parallel
)

if micromamba env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  micromamba install -y -n "$ENV_NAME" -c conda-forge -c bioconda "${PKGS[@]}"
else
  micromamba create  -y -n "$ENV_NAME" -c conda-forge -c bioconda "${PKGS[@]}"
fi

micromamba run -n "$ENV_NAME" bash -lc 'metabat2 -h >/dev/null; snakemake --version >/dev/null'
