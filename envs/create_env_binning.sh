#!/usr/bin/env bash
# Create/update Binny env and clone workflow
set -euo pipefail

# ── Paths
ROOT="/media/Box/MeMaMoRe"
ENV_NAME="binny_env"
BINNY_DIR="$ROOT/external/binny"

# ── Micromamba bootstrap
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
eval "$(micromamba shell hook --shell=bash)"

# ── Create or update env
if micromamba env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  micromamba install -y -n "$ENV_NAME" -c conda-forge -c bioconda \
    python=3.11 snakemake bedtools samtools pigz yq
else
  micromamba create -y -n "$ENV_NAME" -c conda-forge -c bioconda \
    python=3.11 snakemake bedtools samtools pigz yq
fi

# ── Clone Binny if missing
mkdir -p "$ROOT/external"
if [[ ! -d "$BINNY_DIR/.git" ]]; then
  git clone https://github.com/a-h-b/binny.git "$BINNY_DIR"
fi

# ── Quick smoke test
micromamba run -n "$ENV_NAME" snakemake --version >/dev/null
echo "OK: $ENV_NAME ready. Binny at $BINNY_DIR"
