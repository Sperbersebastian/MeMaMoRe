#!/usr/bin/env bash
# File: envs/create_env_binning.sh
set -euo pipefail

# -----------------------------
# Micromamba bootstrap
# -----------------------------
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
command -v micromamba >/dev/null || { echo "[binning] micromamba not found"; exit 127; }
eval "$(micromamba shell hook --shell=bash)"

# -----------------------------
# Paths & constants
# -----------------------------
ROOT="${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
ENV_NAME="env_binning"
COMEBIN_ENV="comebin_env"

BINNY_FORK_URL="https://github.com/Sperbersebastian/binny"
BINNY_DIR="$ROOT/external/binny"
BINNY_ENVS_PREFIX="$BINNY_DIR/conda/envs"

MAGS_DIR="$ROOT/external/MAGScoT"
MANTIS_REF="$ROOT/refdata/mantis/NCBI"

# Core packages for env_binning
PKGS=(
  snakemake conda conda-libmamba-solver mamba
  yq
  metabat2
  checkm2 gtdbtk
  barrnap tRNAscan-SE
  samtools seqkit
  r-base r-optparse r-dplyr r-readr r-funr r-digest
  hmmer prodigal parallel git curl
  python
)

# -----------------------------
# Create or update env_binning
# -----------------------------
if micromamba env list | awk 'NR>2{print $1}' | grep -qx "$ENV_NAME"; then
  echo "[binning] updating $ENV_NAME"
  micromamba install -y -n "$ENV_NAME" -c conda-forge -c bioconda "${PKGS[@]}"
else
  echo "[binning] creating $ENV_NAME"
  micromamba create  -y -n "$ENV_NAME" -c conda-forge -c bioconda "${PKGS[@]}"
fi

# Sanity
micromamba run -n "$ENV_NAME" bash -lc '
  snakemake --version >/dev/null
  metabat2 -h        >/dev/null
  Rscript --version  >/dev/null
'
echo "[ok] core tools present in $ENV_NAME"

# -----------------------------
# COMEBin env (auto CPU/GPU)
# -----------------------------
if micromamba env list | awk 'NR>2{print $1}' | grep -qx "$COMEBIN_ENV"; then
  echo "[comebin] env exists: $COMEBIN_ENV"
else
  echo "[comebin] creating $COMEBIN_ENV (auto GPU detect)"
  if command -v nvidia-smi >/dev/null 2>&1; then
    echo "[gpu] NVIDIA GPU detected -> CUDA build"
    micromamba create -y -n "$COMEBIN_ENV" -c bioconda -c conda-forge -c pytorch -c nvidia \
      comebin pytorch torchvision torchaudio "pytorch-cuda=11.8"
  else
    echo "[cpu] No GPU -> CPU-only build"
    micromamba create -y -n "$COMEBIN_ENV" -c bioconda -c conda-forge -c pytorch \
      comebin pytorch torchvision torchaudio cpuonly
  fi
fi

# PyTorch info
micromamba run -n "$COMEBIN_ENV" python - <<'PY' || true
import torch
print("PyTorch:", torch.__version__, "CUDA available:", torch.cuda.is_available())
PY

# COMEBin CLI (Groß/Klein)
COMEBIN_CLI="$(micromamba run -n "$COMEBIN_ENV" bash -lc 'command -v COMEBin || command -v comebin || true')"
if [[ -n "${COMEBIN_CLI:-}" ]]; then
  echo "[ok] COMEBin CLI at: $COMEBIN_CLI"
  micromamba run -n "$COMEBIN_ENV" "$COMEBIN_CLI" -h | head -n 3 || true
else
  echo "[warn] COMEBin CLI not found in $COMEBIN_ENV"
  micromamba run -n "$COMEBIN_ENV" bash -lc 'ls -1 "$CONDA_PREFIX/bin" || true'
fi

# -----------------------------
# MAGScoT (repo only)
# -----------------------------
if [[ ! -d "$MAGS_DIR/.git" ]]; then
  echo "[magscot] cloning"
  git clone https://github.com/ikmb/MAGScoT "$MAGS_DIR"
else
  echo "[magscot] updating"
  git -C "$MAGS_DIR" fetch --all --prune
  git -C "$MAGS_DIR" pull --ff-only || true
fi
echo "[ok] MAGScoT at: $MAGS_DIR"

# -----------------------------
# Binny fork with detached-HEAD fix
# -----------------------------
mkdir -p "$BINNY_DIR"
if [[ ! -d "$BINNY_DIR/.git" ]]; then
  echo "[binny] cloning fork"
  git clone "$BINNY_FORK_URL" "$BINNY_DIR"
else
  echo "[binny] updating fork"
  git -C "$BINNY_DIR" fetch --all --prune
  # if detached, switch to main tracking origin/main
  if ! git -C "$BINNY_DIR" symbolic-ref -q HEAD >/dev/null; then
    git -C "$BINNY_DIR" checkout -B main origin/main || true
  fi
  git -C "$BINNY_DIR" branch -u origin/main main 2>/dev/null || true
  git -C "$BINNY_DIR" checkout main || true
  git -C "$BINNY_DIR" pull --ff-only || true
fi

# Ensure rule envs prefix exists
mkdir -p "$BINNY_ENVS_PREFIX"

# -----------------------------
# Patch Binny config: envs + in-path snakemake
# -----------------------------
(
  cd "$BINNY_DIR"
  sed -i -e "s|conda_source: \"\"|conda_source: \"${BINNY_ENVS_PREFIX}\"|g" \
         -e "s|snakemake_env: \"\"|snakemake_env: \"in_path\"|g" config/config.*.yaml
  # prokka_env/mantis_env leer lassen -> werden unter conda_source erzeugt
)
echo "[binny] configs patched: conda_source='${BINNY_ENVS_PREFIX}', snakemake_env='in_path'"

# -----------------------------
# MANTIS: NCBI gc.prt Referenz
# -----------------------------
mkdir -p "$MANTIS_REF"
if [[ ! -s "$MANTIS_REF/gc.prt.dmp" ]]; then
  echo "[mantis] fetching NCBI gc.prt"
  curl -fsSL https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt -o "$MANTIS_REF/gc.prt.dmp"
fi
cp -f "$MANTIS_REF/gc.prt.dmp" "$MANTIS_REF/gc.prt"
: > "$MANTIS_REF/README"

CFG_FILE="$BINNY_DIR/config/binny_mantis.cfg"
if [[ -f "$CFG_FILE" ]] && grep -q '^ncbi_ref_folder=' "$CFG_FILE"; then
  sed -i -E "s|^ncbi_ref_folder=.*|ncbi_ref_folder=${MANTIS_REF}/|" "$CFG_FILE"
else
  printf "\nncbi_ref_folder=%s/\n" "$MANTIS_REF" >> "$CFG_FILE"
fi
echo "[mantis] binny_mantis.cfg -> ncbi_ref_folder=$MANTIS_REF/"

# -----------------------------
# Optional: mirror gc.prt in bereits gebaute Rule-Envs
# -----------------------------
if compgen -G "${BINNY_ENVS_PREFIX}/*" >/dev/null; then
  while IFS= read -r -d '' envp; do
    [[ -d "$envp" ]] || continue
    py_site="$(micromamba run -p "$envp" python - <<'PY'
import sysconfig; print(sysconfig.get_paths().get("purelib") or "")
PY
)" || py_site=""
    [[ -n "$py_site" ]] || continue
    dst="$py_site/Resources/NCBI"; mkdir -p "$dst"
    ln -sf "$MANTIS_REF/gc.prt.dmp" "$dst/gc.prt.dmp"
    ln -sf "$MANTIS_REF/gc.prt.dmp" "$dst/gc.prt"
    : > "$py_site/Resources/Taxonomy.db"
    echo "[mantis] linked gc.prt* -> $dst/"
  done < <(find "$BINNY_ENVS_PREFIX" -maxdepth 1 -type d -print0)
else
  echo "[mantis] note: rule envs will be created by Snakemake later."
fi

# -----------------------------
# Print run helper
# -----------------------------
cat <<HELP

[ready] Binning stack prepared.

Run Binny:

  export SNAKEMAKE_CONDA_FRONTEND=mamba
  export MAMBA_EXE="\$(command -v micromamba)"
  export SNAKEMAKE_CONDA_PREFIX="$BINNY_ENVS_PREFIX"
  micromamba run -n "$ENV_NAME" bash -lc "cd \"$BINNY_DIR\" && ./binny -l -n 'TESTRUN' -r config/config.test.yaml"

COMEBin env: $COMEBIN_ENV
MAGScoT path: $MAGS_DIR
HELP
