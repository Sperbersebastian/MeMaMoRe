#!/usr/bin/env bash
# Create/update binning envs and assets: Binny driver + COMEBin
set -Eeuo pipefail

# ── Paths
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BINNY_DIR="$ROOT/external/binny"
COMEBIN_DIR="$ROOT/external/comebin"
SNK_PREFIX="${SNAKEMAKE_CONDA_PREFIX:-$HOME/micromamba/envs/_snk_binny_persist}"

# micromamba + deterministic Snakemake conda prefix
# ── Micromamba bootstrap (force our root + binary)
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MM="$MAMBA_ROOT_PREFIX/bin/micromamba"
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
command -v "$MM" >/dev/null || { echo "micromamba binary not found at $MM"; exit 1; }
eval "$("$MM" shell hook --shell=bash)"

# Fix Snakemake prefix. Avoid inherited bad values like /envs
unset SNAKEMAKE_CONDA_PREFIX
SNK_PREFIX="${SNAKEMAKE_CONDA_PREFIX:-$HOME/micromamba/envs/_snk_binny_persist}"
export SNAKEMAKE_CONDA_PREFIX="$SNK_PREFIX"
mkdir -p "$SNAKEMAKE_CONDA_PREFIX"



############################################
# 1) Binny driver env (adds metabat2)
############################################
ENV_BINNY="binny_env"
if "$MM" env list | awk '{print $1}' | grep -qx "$ENV_BINNY"; then
  "$MM" install -y -n "$ENV_BINNY" -c conda-forge -c bioconda \
    python=3.11 snakemake hmmer diamond mantis_pfa prodigal bedtools samtools pigz yq nltk metabat2 'conda>=24,<25'
else
  "$MM" create  -y -n "$ENV_BINNY" -c conda-forge -c bioconda \
    python=3.11 snakemake hmmer diamond mantis_pfa prodigal bedtools samtools pigz yq nltk metabat2 'conda>=24,<25'
fi
micromamba activate "$ENV_BINNY"

# Fixed prefix for Snakemake sub-envs
export SNAKEMAKE_CONDA_PREFIX="$SNK_PREFIX"

# Binny checkout
if [[ ! -d "$BINNY_DIR" ]]; then
  mkdir -p "$ROOT/external"
  git clone https://github.com/Sperbersebastian/binny.git "$BINNY_DIR"
fi

# Vendor tiny assets for MANTIS
ROOT="$ROOT" bash "$ROOT/modules/binning/binny/install.sh"

# Pre-create MANTIS Snakemake env
pushd "$BINNY_DIR" >/dev/null
snakemake -s Snakefile --cores 1 --use-conda \
  --conda-prefix "$SNK_PREFIX" --conda-create-envs-only -R mantis_checkm_marker_sets
popd >/dev/null

# Inject runtime files + cfg wiring
MANTIS_BIN="$(find "$SNK_PREFIX" -type f -path '*/bin/mantis' | head -n1)"
[[ -x "$MANTIS_BIN" ]] || { echo "mantis not found under $SNK_PREFIX"; exit 2; }
MENV="$(readlink -f "$(dirname "$MANTIS_BIN")/..")"
PYVER="$(basename "$(find "$MENV/lib" -maxdepth 1 -type d -name 'python3*' | head -n1)")"
DEST_RES="$MENV/lib/$PYVER/site-packages/Resources"
install -D -m0644 "$BINNY_DIR/database/support/Resources/NCBI/gc.prt.dmp" "$DEST_RES/NCBI/gc.prt.dmp"
sed -i -E \
  -e "s|^custom_ref=.*/checkm_tf/checkm_filtered_tf\.hmm$|custom_ref=$BINNY_DIR/database/hmms/checkm_tf/checkm_filtered_tf.hmm|" \
  -e "s|^custom_ref=.*/checkm_pf/checkm_filtered_pf\.hmm$|custom_ref=$BINNY_DIR/database/hmms/checkm_pf/checkm_filtered_pf.hmm|" \
  "$BINNY_DIR/config/binny_mantis.cfg"
export NLTK_DATA="$BINNY_DIR/database/support/nltk_data"
"$MANTIS_BIN" check_sql -mc "$BINNY_DIR/config/binny_mantis.cfg"
"$MANTIS_BIN" check     -mc "$BINNY_DIR/config/binny_mantis.cfg" --no_taxonomy || true

echo "OK: $ENV_BINNY ready. Binny at $BINNY_DIR. Snakemake env prefix: $SNK_PREFIX"

############################################
# 2) COMEBin env (Bioconda, Py 3.7)
############################################
ENV_COMEBIN="comebin_env"
if "$MM" env list | awk '{print $1}' | grep -qx "$ENV_COMEBIN"; then
  "$MM" install -y -n "$ENV_COMEBIN" --channel-priority flexible \
    -c bioconda -c conda-forge -c defaults \
    python=3.7 "comebin>=1.0.4"
else
  "$MM" create  -y -n "$ENV_COMEBIN" --channel-priority flexible \
    -c bioconda -c conda-forge -c defaults \
    python=3.7 "comebin>=1.0.4"
fi
echo "OK: base $ENV_COMEBIN created/updated from Bioconda (Py 3.7)."

# Fallback repo checkout for runner script if bioconda lacks it
mkdir -p "$ROOT/external"
COMEBIN_DIR="$ROOT/external/comebin"
if [[ ! -d "$COMEBIN_DIR/.git" ]]; then
  git clone https://github.com/ziyewang/COMEBin.git "$COMEBIN_DIR"
fi
if "$MM" run -n "$ENV_COMEBIN" bash -lc 'command -v run_comebin.sh >/dev/null || command -v comebin >/dev/null'; then
  echo "OK: COMEBin runner available in $ENV_COMEBIN."
else
  echo "Note: runner not on PATH in env; pipeline will use $COMEBIN_DIR/run_comebin.sh."
fi


echo "OK: $ENV_COMEBIN ready."
