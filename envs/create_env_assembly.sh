#!/usr/bin/env bash
set -euo pipefail
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
eval "$(micromamba shell hook --shell=bash)"


CH="-c conda-forge -c bioconda --channel-priority flexible"

# env_assembly_core → SPAdes (and python)
if micromamba env list | awk 'NR>2{print $1}' | grep -qx env_assembly_core; then
  micromamba install -y -n env_assembly_core $CH python=3.10 spades
else
  micromamba create  -y -n env_assembly_core $CH python=3.10 spades
fi

# env_mapping_coverm → mapping + coverage tools
if micromamba env list | awk 'NR>2{print $1}' | grep -qx env_mapping_coverm; then
  micromamba install -y -n env_mapping_coverm $CH python=3.10 bwa-mem2 samtools coverm
else
  micromamba create  -y -n env_mapping_coverm $CH python=3.10 bwa-mem2 samtools coverm
fi

echo "[creator] ensured: env_assembly_core, env_mapping_coverm"
