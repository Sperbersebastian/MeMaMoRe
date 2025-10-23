#!/usr/bin/env bash
set -euo pipefail
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
eval "$(micromamba shell hook --shell=bash)"

#!/usr/bin/env bash
set -euo pipefail

# QC core: trimming + QC
micromamba create -y -n env_qc_core \
  -c bioconda -c conda-forge \
  --channel-priority flexible \
  python=3.10 fastp fastqc

# QC aggregation: MultiQC
micromamba create -y -n env_qc_mqc \
  -c bioconda -c conda-forge \
  --channel-priority flexible \
  python=3.9 multiqc=1.17
