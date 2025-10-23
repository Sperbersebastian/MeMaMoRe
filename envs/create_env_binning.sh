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
ROOT="${ROOT:-/media/Box/MeMaMoRe}"
ENV_NAME="env_binning"

BINNY_FORK_URL="https://github.com/Sperbersebastian/binny"
BINNY_DIR="$ROOT/external/binny"
BINNY_ENVS_PREFIX="$BINNY_DIR/conda/envs"   # Snakemake rule envs live here

COMEBIN_ENV="comebin_env"
MAGS_DIR="$ROOT/external/MAGScoT"

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

# Quick sanity
micromamba run -n "$ENV_NAME" bash -lc '
  snakemake --version >/dev/null
  metabat2 -h        >/dev/null
  Rscript --version  >/dev/null
'
echo "[ok] core tools present in $ENV_NAME"

# -----------------------------
# COMEBin in its own env
# -----------------------------
if micromamba env list | awk 'NR>2{print $1}' | grep -qx "$COMEBIN_ENV"; then
  echo "[comebin] env exists: $COMEBIN_ENV"
else
  echo "[comebin] creating $COMEBIN_ENV"
  micromamba create -y -n "$COMEBIN_ENV" -c conda-forge -c bioconda -c pytorch -c nvidia \
    comebin "pytorch-cuda=11.8" pytorch
fi

# test GPU note (optional)
micromamba run -n "$COMEBIN_ENV" python - <<'PY' || true
import torch
print(f"PyTorch: {torch.__version__} CUDA available:", torch.cuda.is_available())
PY

# Try COMEBin CLI (if wrapper present)
if micromamba run -n "$COMEBIN_ENV" bash -lc 'command -v comebin >/dev/null'; then
  echo "[ok] COMEBin CLI:"
  micromamba run -n "$COMEBIN_ENV" comebin -h | head -n 3 || true
else
  echo "[warn] comebin CLI not found in $COMEBIN_ENV"
fi

# -----------------------------
# MAGScoT (no env; uses env_binning tools)
# -----------------------------
if [[ ! -d "$MAGS_DIR" ]]; then
  echo "[magscot] cloning"
  git clone https://github.com/ikmb/MAGScoT "$MAGS_DIR"
else
  echo "[magscot] updating"
  git -C "$MAGS_DIR" fetch --all --prune
  git -C "$MAGS_DIR" pull --ff-only || true
fi
echo "[ok] MAGScoT at: $MAGS_DIR"

# -----------------------------
# Binny fork
# -----------------------------
mkdir -p "$BINNY_DIR"
if [[ ! -d "$BINNY_DIR/.git" ]]; then
  echo "[binny] cloning fork"
  git clone "$BINNY_FORK_URL" "$BINNY_DIR"
else
  echo "[binny] updating fork"
  git -C "$BINNY_DIR" fetch --all --prune
  if ! git -C "$BINNY_DIR" rev-parse --abbrev-ref --symbolic-full-name '@{u}' >/dev/null 2>&1; then
    if git -C "$BINNY_DIR" rev-parse --verify origin/main >/dev/null 2>&1; then
      git -C "$BINNY_DIR" branch --set-upstream-to=origin/main || true
    fi
  fi
  git -C "$BINNY_DIR" pull --rebase || true
fi

# Ensure rule envs prefix exists
mkdir -p "$BINNY_ENVS_PREFIX"

# -----------------------------
# Patch Binny config to use our envs & in-path Snakemake
# -----------------------------
(
  cd "$BINNY_DIR"
  my_conda_env_path="$BINNY_ENVS_PREFIX"
  my_prokka_env=""      # let Binny manage
  my_mantis_env=""      # let Binny manage
  my_snakemake_env="in_path"

  sed -i -e "s|conda_source: \"\"|conda_source: \"${my_conda_env_path}\"|g" config/config.*.yaml
  sed -i -e "s|prokka_env: \"\"|prokka_env: \"${my_prokka_env}\"|g" \
         -e "s|mantis_env: \"\"|mantis_env: \"${my_mantis_env}\"|g" \
         -e "s|snakemake_env: \"\"|snakemake_env: \"${my_snakemake_env}\"|g" config/config.*.yaml
)
echo "[binny] configs patched: conda_source='${BINNY_ENVS_PREFIX}', snakemake_env='in_path'"

# -----------------------------
# MANTIS: NCBI gc.prt.dmp + CFG update
# -----------------------------
MANTIS_REF="$ROOT/refdata/mantis/NCBI"
mkdir -p "$MANTIS_REF"
if [[ ! -s "$MANTIS_REF/gc.prt.dmp" ]]; then
  echo "[mantis] fetching NCBI gc.prt"
  curl -fsSL https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt -o "$MANTIS_REF/gc.prt.dmp"
fi
# zweiter Name, damit der Checker beides findet
cp -f "$MANTIS_REF/gc.prt.dmp" "$MANTIS_REF/gc.prt"
: > "$MANTIS_REF/README"

# CFG so setzen, dass Snakefile nur ausliest, aber nichts hart codiert
CFG_FILE="$BINNY_DIR/config/binny_mantis.cfg"
if grep -q '^ncbi_ref_folder=' "$CFG_FILE"; then
  sed -i -E "s|^ncbi_ref_folder=.*|ncbi_ref_folder=${MANTIS_REF}/|" "$CFG_FILE"
else
  printf "\nncbi_ref_folder=%s/\n" "$MANTIS_REF" >> "$CFG_FILE"
fi
echo "[mantis] binny_mantis.cfg -> ncbi_ref_folder=$MANTIS_REF/"


# -----------------------------
# OPTIONAL: NLTK tagger seed (for mantis envs)
# -----------------------------
seed_nltk() {
  local env_path="$1"
  echo "[mantis] seeding NLTK tagger in: $env_path"
  micromamba run -p "$env_path" bash -lc '
    set -e
    export NLTK_DATA="$CONDA_PREFIX/share/nltk_data"
    python - <<PY
import nltk, os
d=os.environ.get("NLTK_DATA","")
os.makedirs(d, exist_ok=True)
try:
    nltk.data.find("taggers/averaged_perceptron_tagger_eng")
except LookupError:
    nltk.download("averaged_perceptron_tagger_eng", download_dir=d)
print("NLTK_DATA =", d)
PY
  ' || true
}

# -----------------------------
# OPTION A: mirror gc.prt.dmp into any existing rule envs
# -----------------------------
mirror_gc_to_env() {
  local env_path="$1"
  local py_site
  py_site="$(micromamba run -p "$env_path" python - <<'PY'
import sysconfig
print(sysconfig.get_paths().get("purelib") or "")
PY
)" || py_site=""
  if [[ -z "$py_site" ]]; then
    echo "[mantis] WARN: cannot locate site-packages for $env_path"
    return 0
  fi
  local dst_dir="$py_site/Resources/NCBI"
  mkdir -p "$dst_dir"
  ln -sf "$MANTIS_REF/gc.prt.dmp" "$dst_dir/gc.prt.dmp"
  ln -sf "$MANTIS_REF/gc.prt.dmp" "$dst_dir/gc.prt"
  # lightweight presence file
  : > "${py_site}/Resources/Taxonomy.db"
  echo "[mantis] linked gc.prt(.dmp) -> $dst_dir/"
}

# Run seeding on any already-built rule envs
if compgen -G "${BINNY_ENVS_PREFIX}/*" >/dev/null; then
  while IFS= read -r -d '' envp; do
    [[ -d "$envp" ]] || continue
    [[ "$envp" =~ mantis ]] && seed_nltk "$envp" || true
    mirror_gc_to_env "$envp"
  done < <(find "$BINNY_ENVS_PREFIX" -maxdepth 1 -type d -print0)
else
  echo "[mantis] note: rule envs will be created by Snakemake; re-run this script after first build to mirror gc.prt if nötig. Alternativ übernimmt das Snakefile den Pre-Run-Seeding-Schritt."
fi

# -----------------------------
# Print run helper
# -----------------------------
cat <<HELP

[ready] Binning stack prepared.

Run Binny using your shared env_binning:

  export SNAKEMAKE_CONDA_FRONTEND=mamba
  export MAMBA_EXE="$(command -v micromamba)"
  export SNAKEMAKE_CONDA_PREFIX="$BINNY_ENVS_PREFIX"

  micromamba run -n "$ENV_NAME" bash -lc "cd \"$BINNY_DIR\" && ./binny -l -n 'TESTRUN' -r config/config.test.yaml"

If MANTIS rule envs are built later, re-run this script to mirror gc.prt into their site-packages/Resources/NCBI.
Snakefile kann zusätzlich vor jedem 'mantis run' seeden, ohne feste Pfade.

COMEBin environment: $COMEBIN_ENV
MAGScoT path:       $MAGS_DIR

HELP
