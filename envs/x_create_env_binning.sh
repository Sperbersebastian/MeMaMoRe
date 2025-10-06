#!/usr/bin/env bash
# Bootstrap & run: MetaBAT2 → Binny → (optional) COMEBin/MAGScoT → CheckM2 → GTDB-Tk → rRNA/tRNA
# Usage: ./run.sh CONTIGS.fa [BAM ...]     # BAMs must be coordinate-sorted and indexed (.bai)
set -euo pipefail
[[ "${DEBUG:-0}" == "1" ]] && set -x

# ── Repo + micromamba ────────────────────────────────────────────────────────
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"

if command -v micromamba >/dev/null 2>&1; then MM="$(command -v micromamba)"
elif [[ -x "$MAMBA_ROOT_PREFIX/bin/micromamba" ]]; then MM="$MAMBA_ROOT_PREFIX/bin/micromamba"
else echo "micromamba not found. Install it or set MAMBA_ROOT_PREFIX." >&2; exit 127; fi

# Snakemake uses micromamba for rule envs
export SNAKEMAKE_CONDA_FRONTEND=micromamba
export SNAKEMAKE_CONDA_PREFIX="${SNAKEMAKE_CONDA_PREFIX:-$REPO_ROOT/.snk_envs/_snk_binny_persist}"
mkdir -p "$SNAKEMAKE_CONDA_PREFIX"

# ── Inputs ───────────────────────────────────────────────────────────────────
if [[ $# -lt 1 ]]; then echo "Usage: $(basename "$0") CONTIGS.fa [BAM ...]" >&2; exit 2; fi
CONTIGS="$(readlink -f -- "$1")"; shift || true
BAMS=( "$@" )  # 0..n BAMs
[[ -s "$CONTIGS" ]] || { echo "Missing contigs: $CONTIGS" >&2; exit 2; }

CPUS="${CPUS:-8}"
OUTDIR="${OUTDIR:-$(pwd)/run}"
SAMPLE="${SAMPLE:-sample1}"
EXTERNAL_DEPTH="${CONTIG_DEPTH:-}"   # optional depth table path

# ── Paths ────────────────────────────────────────────────────────────────────
BINNY_DIR="$(cd "$REPO_ROOT/external/binny" && pwd)"
GTDBTK_DATA_PATH="${GTDBTK_DATA_PATH:-$REPO_ROOT/refdata/gtdbtk/release226}"
[[ -d "$GTDBTK_DATA_PATH" ]] && export GTDBTK_DATA_PATH || true

mkdir -p \
  "$OUTDIR"/{depth,metabat2,binny,comebin,magscot,checkm2,gtdbtk,ann/barrnap,ann/trnascan} \
  "$REPO_ROOT/logs"

mm_run(){ "$MM" run -n "$ENV_MAIN" "$@"; }
env_has(){ "$MM" run -n "$ENV_MAIN" bash -lc "command -v '$1' >/dev/null"; }

# ── Env bootstrap ────────────────────────────────────────────────────────────
ENV_MAIN="${ENV_MAIN:-env_binning}"
REQ_PKGS=(metabat2 snakemake yq checkm2 gtdbtk barrnap tRNAscan-SE samtools)

env_exists() { "$MM" env list | awk '{print $1}' | grep -qx "$ENV_MAIN"; }
have_in_env(){ "$MM" run -n "$ENV_MAIN" bash -lc "command -v '$1' >/dev/null"; }
pkg_bin(){ case "$1" in metabat2) echo jgi_summarize_bam_contig_depths;; tRNAscan-SE) echo tRNAscan-SE;; *) echo "$1";; esac; }

ensure_env(){
  if ! env_exists; then
    echo "[bootstrap] creating $ENV_MAIN"
    "$MM" create -y -n "$ENV_MAIN" -c conda-forge -c bioconda "${REQ_PKGS[@]}"
  else
    missing=(); for p in "${REQ_PKGS[@]}"; do b="$(pkg_bin "$p")"; have_in_env "$b" || missing+=("$p"); done
    ((${#missing[@]})) && "$MM" install -y -n "$ENV_MAIN" -c conda-forge -c bioconda "${missing[@]}"
  fi
}
ensure_env

# ── 0) Coverage table ────────────────────────────────────────────────────────
DEPTH_FILE=""
if (( ${#BAMS[@]} > 0 )); then
  for b in "${BAMS[@]}"; do
    [[ -s "$b" ]] || { echo "BAM missing: $b" >&2; exit 2; }
    [[ -s "${b}.bai" || -s "${b%.bam}.bai" ]] || { echo "Index (.bai) missing for: $b" >&2; exit 2; }
  done
  mm_run jgi_summarize_bam_contig_depths --outputDepth "$OUTDIR/depth/depth.txt" "${BAMS[@]}"
  [[ -s "$OUTDIR/depth/depth.txt" ]] || { echo "Depth file empty" >&2; exit 2; }
  DEPTH_FILE="$OUTDIR/depth/depth.txt"
elif [[ -n "$EXTERNAL_DEPTH" ]]; then
  DEPTH_FILE="$(readlink -f -- "$EXTERNAL_DEPTH")"
  [[ -s "$DEPTH_FILE" ]] || { echo "Provided CONTIG_DEPTH not found: $DEPTH_FILE" >&2; exit 2; }
else
  echo "No BAMs or CONTIG_DEPTH provided. MetaBAT2 can run without depth, Binny cannot." >&2
fi

# ── 1) MetaBAT2 ──────────────────────────────────────────────────────────────
DEPTH_ARG=(); [[ -n "$DEPTH_FILE" ]] && DEPTH_ARG=( -a "$DEPTH_FILE" )
mm_run metabat2 -i "$CONTIGS" "${DEPTH_ARG[@]}" -t "$CPUS" -o "$OUTDIR/metabat2/${SAMPLE}.bin"

# Pick bins dir for downstream (MetaBAT2 output pattern)
BINS_DIR="$OUTDIR/metabat2"

# ── 2) Binny (Snakemake with micromamba rule envs) ───────────────────────────
BINNY_CFG="$OUTDIR/binny/config.yaml"; mkdir -p "$(dirname "$BINNY_CFG")"

if [[ -n "$DEPTH_FILE" ]]; then
  cat > "$BINNY_CFG" <<EOF
raws:
  assembly: "$CONTIGS"
  metagenomics_alignment: ""
  contig_depth: "$DEPTH_FILE"
sample: "$SAMPLE"
outputdir: "$OUTDIR/binny"
conda_source: "$SNAKEMAKE_CONDA_PREFIX"
EOF
elif (( ${#BAMS[@]} > 0 )); then
  first_bam="$(readlink -f -- "${BAMS[0]}")"
  cat > "$BINNY_CFG" <<EOF
raws:
  assembly: "$CONTIGS"
  metagenomics_alignment: "$first_bam"
  contig_depth: ""
sample: "$SAMPLE"
outputdir: "$OUTDIR/binny"
conda_source: "$SNAKEMAKE_CONDA_PREFIX"
EOF
else
  echo "Binny requires BAM or depth. Skipping Binny." >&2
fi

if [[ -f "$BINNY_CFG" ]]; then
  mm_run snakemake -s "$BINNY_DIR/Snakefile" \
    --cores "$CPUS" \
    --use-conda \
    --conda-prefix "$SNAKEMAKE_CONDA_PREFIX" \
    --conda-frontend "$SNAKEMAKE_CONDA_FRONTEND" \
    --configfile "$BINNY_CFG" -p all
fi

# ── 3) COMEBin (optional) ────────────────────────────────────────────────────
if env_has run_comebin.sh; then
  if (( ${#BAMS[@]} > 0 )); then
    mm_run run_comebin.sh -r "$CONTIGS" -b "${BAMS[@]}" -o "$OUTDIR/comebin" -t "$CPUS"
  else
    mm_run run_comebin.sh -r "$CONTIGS" -o "$OUTDIR/comebin" -t "$CPUS"
  fi
elif env_has comebin; then
  if (( ${#BAMS[@]} > 0 )); then
    mm_run comebin -r "$CONTIGS" -b "${BAMS[@]}" -o "$OUTDIR/comebin" -t "$CPUS"
  else
    mm_run comebin -r "$CONTIGS" -o "$OUTDIR/comebin" -t "$CPUS"
  fi
else
  echo "COMEBin not installed — skipping."
fi

# ── 4) MAGScoT (optional) ────────────────────────────────────────────────────
if env_has magscot || env_has MAGScot; then
  mm_run magscot -i "$CONTIGS" -o "$OUTDIR/magscot" -t "$CPUS" || true
else
  echo "MAGScoT not installed — skipping."
fi

# ── 5) CheckM2 ───────────────────────────────────────────────────────────────
if env_has checkm2; then
  mm_run checkm2 predict --input "$BINS_DIR" --extension fa -t "$CPUS" --output-directory "$OUTDIR/checkm2"
else
  echo "CheckM2 not installed — skipping."
fi

# ── 6) GTDB-Tk (only if data path exists) ────────────────────────────────────
if env_has gtdbtk && [[ -n "${GTDBTK_DATA_PATH:-}" && -d "$GTDBTK_DATA_PATH" ]]; then
  mm_run gtdbtk classify_wf --genome_dir "$BINS_DIR" --out_dir "$OUTDIR/gtdbtk" --cpus "$CPUS" --skip_ani_screen
else
  echo "GTDB-Tk or GTDBTK_DATA_PATH missing — skipping."
fi

# ── 7/8) Barrnap + tRNAscan-SE ───────────────────────────────────────────────
find "$BINS_DIR" -maxdepth 1 -type f \( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' \) |
while read -r f; do
  base="$(basename "${f%.*}")"
  env_has barrnap && mm_run barrnap --threads "$CPUS" < "$f" > "$OUTDIR/ann/barrnap/${base}.gff" || true
  env_has tRNAscan-SE && mm_run tRNAscan-SE -o "$OUTDIR/ann/trnascan/${base}.txt" -f "$OUTDIR/ann/trnascan/${base}.ss" "$f" || true
done

echo "DONE. Outputs under: $OUTDIR"
