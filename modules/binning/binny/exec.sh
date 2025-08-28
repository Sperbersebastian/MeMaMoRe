#!/usr/bin/env bash
# Robust Binny runner (self-bootstraps MANTIS assets & NLTK)
# expects env: ROOT PARAMS_YAML SAMPLE CONTIGS BAM
set -Eeuo pipefail
umask 002

# -------- helpers --------
fail(){ echo "[binny] ERROR: $*" >&2; exit 2; }
on_err(){ echo "[binny] crashed at line $1"; }
trap 'on_err $LINENO' ERR

need(){ command -v "$1" >/dev/null || fail "missing executable: $1"; }
abspath(){ readlink -f -- "$1"; }

# -------- inputs --------
[[ -n "${ROOT:-}"        ]] || fail "ROOT not set"
[[ -n "${PARAMS_YAML:-}" ]] || fail "PARAMS_YAML not set"
[[ -n "${SAMPLE:-}"      ]] || fail "SAMPLE not set"
[[ -n "${CONTIGS:-}"     ]] || fail "CONTIGS not set"
[[ -n "${BAM:-}"         ]] || fail "BAM not set"
[[ -s "$CONTIGS" ]] || fail "$SAMPLE no CONTIGS at $CONTIGS"
[[ -s "$BAM"     ]] || fail "$SAMPLE no BAM at $BAM"

ROOT="$(abspath "$ROOT")"
BINNY_DIR="$ROOT/external/binny"
[[ -d "$BINNY_DIR" ]] || fail "missing $BINNY_DIR"

OUTDIR="$ROOT/SRA/binning/binny/$SAMPLE"
LOGS="$ROOT/logs"
mkdir -p "$OUTDIR" "$LOGS"

# -------- micromamba / snakemake driver env --------
# micromamba + deterministic Snakemake conda prefix
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
command -v micromamba >/dev/null || { echo "micromamba missing"; exit 1; }
eval "$(micromamba shell hook --shell=bash)"

# discover root prefix and set a stable sub-env prefix for Snakemake
ROOT_PREFIX="$(micromamba info | awk -F': ' '/Root prefix/ {print $2}')"
_default_snk="$ROOT_PREFIX/envs/_snk_binny_persist"
export SNAKEMAKE_CONDA_PREFIX="${SNAKEMAKE_CONDA_PREFIX:-$_default_snk}"


if ! micromamba env list | grep -qE '^\s*binny_env\s'; then
  micromamba create -y -n binny_env -c conda-forge -c bioconda \
    python=3.11 snakemake hmmer diamond mantis_pfa prodigal "nltk>=3.8" \
    conda "conda-libmamba-solver>=25"
fi
micromamba activate binny_env
need snakemake
need conda
conda info --json >/dev/null || fail "'conda info' failed; reinstall driver env"

# Snakemake sub-env prefix (persistent)
export SNAKEMAKE_CONDA_PREFIX="${SNAKEMAKE_CONDA_PREFIX:-$HOME/micromamba/envs/_snk_binny_persist}"

# -------- ensure vendored assets exist (idempotent) --------
if [[ ! -s "$BINNY_DIR/database/support/Resources/NCBI/gc.prt.dmp" ]] \
 || [[ ! -d "$BINNY_DIR/database/support/nltk_data" ]] \
 || [[ ! -s "$BINNY_DIR/database/hmms/checkm_tf/metadata.tsv" ]] \
 || [[ ! -s "$BINNY_DIR/database/hmms/checkm_pf/metadata.tsv" ]]; then
  [[ -x "$ROOT/modules/binning/binny/install.sh" ]] || fail "installer missing: $ROOT/modules/binning/binny/install.sh"
  ROOT="$ROOT" bash "$ROOT/modules/binning/binny/install.sh"
fi

# -------- pre-create MANTIS env + inject tiny runtime files --------
bootstrap_mantis_env() {
  local snk_prefix="$SNAKEMAKE_CONDA_PREFIX"
  if ! find "$snk_prefix" -type f -path '*/bin/mantis' | grep -q .; then
    pushd "$BINNY_DIR" >/dev/null
    snakemake -s Snakefile --cores 1 --use-conda \
      --conda-prefix "$snk_prefix" --conda-create-envs-only -R mantis_checkm_marker_sets
    popd >/dev/null
  fi

  local mantis_bin mantis_env pyver res_ncbi cfg
  mantis_bin="$(find "$snk_prefix" -type f -path '*/bin/mantis' | head -n1)" || true
  [[ -x "$mantis_bin" ]] || fail "cannot find MANTIS in $snk_prefix"
  mantis_env="$(readlink -f "$(dirname "$mantis_bin")/..")"
  pyver="$(basename "$(find "$mantis_env/lib" -maxdepth 1 -type d -name 'python3*' | head -n1)")"
  res_ncbi="$mantis_env/lib/$pyver/site-packages/Resources/NCBI"
  mkdir -p "$res_ncbi"
  install -m0644 "$BINNY_DIR/database/support/Resources/NCBI/gc.prt.dmp" "$res_ncbi/gc.prt.dmp"
  ln -sf gc.prt.dmp "$res_ncbi/gc.prt"  # harmless safety symlink

  # point NLTK to vendored data (no downloads at runtime)
  export NLTK_DATA="$BINNY_DIR/database/support/nltk_data"

  # ensure config points to THIS checkout
  cfg="$BINNY_DIR/config/binny_mantis.cfg"
  sed -i -E \
    -e "s|^custom_ref=.*/checkm_tf/checkm_filtered_tf\.hmm$|custom_ref=$BINNY_DIR/database/hmms/checkm_tf/checkm_filtered_tf.hmm|" \
    -e "s|^custom_ref=.*/checkm_pf/checkm_filtered_pf\.hmm$|custom_ref=$BINNY_DIR/database/hmms/checkm_pf/checkm_filtered_pf.hmm|" \
    "$cfg"

  # preflight (quiet)
  "$mantis_bin" check_sql -mc "$cfg" >/dev/null
  "$mantis_bin" check     -mc "$cfg" >/dev/null || true
}
bootstrap_mantis_env

# -------- parameters --------
getp(){ micromamba run -n env_sra_tools python "$ROOT/scripts/merge_params.py" --defaults "$PARAMS_YAML" --get "$1"; }
THREADS="$(getp global.threads || echo 8)"
MINLEN="$(getp modules.binning.min_contig_len || echo 1500)"

# -------- per-sample config --------
CFG="$OUTDIR/config.yaml"
cat >"$CFG" <<EOF
samples: ["$SAMPLE"]
assemblies:
  $SAMPLE: "$CONTIGS"
mappings:
  $SAMPLE:
    - "$BAM"
output_dir: "$OUTDIR"
min_contig_len: $MINLEN
threads: $THREADS
EOF

# -------- run binny --------
pushd "$BINNY_DIR" >/dev/null
export NLTK_DATA="${NLTK_DATA:-$BINNY_DIR/database/support/nltk_data}"
export SNAKEMAKE_CONDA_PREFIX
./binny -i "$CFG" &> "$LOGS/binning_binny_${SAMPLE}.log"
popd >/dev/null

# -------- contig→bin map for MAGScoT --------
BINS="$OUTDIR/bins"
MAP_RAW="$OUTDIR/contigs_to_bin.tsv"
MAP_SET="$OUTDIR/contigs_to_bin.with_set.tsv"
: > "$MAP_RAW"
shopt -s nullglob
for fa in "$BINS"/*.fa "$BINS"/*.fasta; do
  bin="$(basename "${fa%.*}")"
  awk -v b="$bin" 'BEGIN{RS=">";FS="\n"} NR>1{split($1,h," "); print b "\t" h[1]}' "$fa" >> "$MAP_RAW"
done
awk '{print $1"\t"$2"\tbinny"}' "$MAP_RAW" > "$MAP_SET"

echo "[binny] $SAMPLE done → $OUTDIR"
