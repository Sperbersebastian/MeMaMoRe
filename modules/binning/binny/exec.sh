#!/usr/bin/env bash
# Run Binny Snakefile via Snakemake + micromamba
# env: ROOT PARAMS_YAML SAMPLE CONTIGS [BAM]
set -euo pipefail

fail(){ echo "[binny] ERROR: $*" >&2; exit 2; }
abspath(){ readlink -f -- "$1"; }

# ---- inputs
[[ -n "${ROOT:-}"        ]] || fail "ROOT not set"
[[ -n "${PARAMS_YAML:-}" ]] || fail "PARAMS_YAML not set"
[[ -n "${SAMPLE:-}"      ]] || fail "SAMPLE not set"
[[ -n "${CONTIGS:-}"     ]] || fail "CONTIGS not set"
[[ -s "$CONTIGS" ]] || fail "$SAMPLE no CONTIGS at $CONTIGS"

BINNY_DIR="$ROOT/external/binny"; [[ -d "$BINNY_DIR" ]] || fail "missing $BINNY_DIR"
OUTDIR="$ROOT/SRA/binning/binny/$SAMPLE"
LOGS="$ROOT/logs"
mkdir -p "$OUTDIR" "$LOGS"

# ---- persistent rule-env prefix (writable)
PERSIST="$ROOT/external/binny/conda/envs"
mkdir -p "$PERSIST"
export SNAKEMAKE_CONDA_PREFIX="$PERSIST"
export SNAKEMAKE_CONDA_FRONTEND=mamba
export MAMBA_EXE="$(command -v micromamba)"

# ---- params (robust defaults if PARAMS_YAML is empty or /dev/null)
CPUS="$(micromamba run -n env_binning yq -r '.modules.binning.cpus // empty' "$PARAMS_YAML" 2>/dev/null || true)"
[[ -n "${CPUS:-}" && "$CPUS" != "null" ]] || CPUS=8

MINLEN="$(micromamba run -n env_binning yq -r '.modules.binning.min_contig_len // empty' "$PARAMS_YAML" 2>/dev/null || true)"
[[ -n "${MINLEN:-}" && "$MINLEN" != "null" ]] || MINLEN=1500

# ---- prefer precomputed depth (strip header)
CDEPTH="$ROOT/SRA/binning/coverage/$SAMPLE/$SAMPLE.depth.txt"
CDEPTH_NOHDR="$ROOT/SRA/binning/coverage/$SAMPLE/$SAMPLE.depth.nohdr.txt"
USE_DEPTH=0
if [[ -s "$CDEPTH" ]]; then
  [[ -s "$CDEPTH_NOHDR" ]] || tail -n +2 "$CDEPTH" > "$CDEPTH_NOHDR"
  USE_DEPTH=1
fi

# ---- MANTIS local assets only (no downloads)
export NLTK_DATA="$BINNY_DIR/database/support/nltk_data"
CFG_MANTIS="$BINNY_DIR/config/binny_mantis.cfg"
sed -i -E \
  -e "s|^custom_ref=.*/checkm_tf/checkm_filtered_tf\.hmm$|custom_ref=$BINNY_DIR/database/hmms/checkm_tf/checkm_filtered_tf.hmm|" \
  -e "s|^custom_ref=.*/checkm_pf/checkm_filtered_pf\.hmm$|custom_ref=$BINNY_DIR/database/hmms/checkm_pf/checkm_filtered_pf.hmm|" \
  "$CFG_MANTIS" || true
grep -q '^checkm_filtered_tf_weight=' "$CFG_MANTIS" || echo "checkm_filtered_tf_weight=0.5" >>"$CFG_MANTIS"
grep -q '^checkm_filtered_pf_weight=' "$CFG_MANTIS" || echo "checkm_filtered_pf_weight=1"   >>"$CFG_MANTIS"
echo "[binny] MANTIS custom-only mode."

# ---- Binny run config
CFG="$OUTDIR/config.yaml"
if [[ "$USE_DEPTH" -eq 1 ]]; then
  cat > "$CFG" <<EOF
raws:
  assembly: "$(abspath "$CONTIGS")"
  metagenomics_alignment: ""
  contig_depth: "$(abspath "$CDEPTH_NOHDR")"
sample: "$SAMPLE"
outputdir: "$(abspath "$OUTDIR")"
conda_source: "$PERSIST"
min_cont_length_cutoff: ${MINLEN}
threads: ${CPUS}
EOF
else
  [[ -n "${BAM:-}" && -s "$BAM" ]] || fail "$SAMPLE needs BAM or depth"
  cat > "$CFG" <<EOF
raws:
  assembly: "$(abspath "$CONTIGS")"
  metagenomics_alignment: "$(abspath "$BAM")"
  contig_depth: ""
sample: "$SAMPLE"
outputdir: "$(abspath "$OUTDIR")"
conda_source: "$PERSIST"
min_cont_length_cutoff: ${MINLEN}
threads: ${CPUS}
EOF
fi

# ---- ensure MANTIS HMMs are pressed and have metadata
prep_hmm_dir() {
  local D="$1"
  local HMM=("$D"/checkm_filtered_*.hmm)
  [[ -s "${HMM[0]}" ]] || fail "[binny] missing HMM in $D (need checkm_filtered_*.hmm)"
  # build indices if any missing
  if ! ls "$D"/checkm_filtered_*.hmm.h3{f,i,m,p} >/dev/null 2>&1; then
    micromamba run -n env_binning hmmpress "${HMM[0]}"
  fi
  # metadata.tsv for MANTIS
  if [[ ! -s "$D/metadata.tsv" ]]; then
    awk -v OFS='\t' '
      /^NAME[ \t]/{n=$2}
      /^ACC[ \t]/{a=$2}
      /^LENG[ \t]/{len=$2; db=(a ~ /^PF/ ? "pfam" : (a ~ /^TIGR/ ? "tigrfam" : "custom"));
                   if(n!=""){print n, "|", "acc:"a, "db:"db, "length:"len; n=a=len=""}}
    ' "${HMM[0]}" > "$D/metadata.tsv"
  fi
}
prep_hmm_dir "$BINNY_DIR/database/hmms/checkm_tf"
prep_hmm_dir "$BINNY_DIR/database/hmms/checkm_pf"


# ---- helper paths for Snakefile scripts
export PYTHONPATH="$BINNY_DIR:$BINNY_DIR/workflow/scripts:${PYTHONPATH:-}"
mkdir -p "$BINNY_DIR/logs"

# ---- run Snakemake inside env_binning
LOG="$LOGS/binning_binny_${SAMPLE}.log"
{
  micromamba run -n env_binning snakemake \
    -s "$BINNY_DIR/Snakefile" \
    --directory "$BINNY_DIR" \
    --cores "$CPUS" \
    --use-conda \
    --conda-prefix "$SNAKEMAKE_CONDA_PREFIX" \
    --conda-frontend "$SNAKEMAKE_CONDA_FRONTEND" \
    --rerun-incomplete \
    --restart-times 1 \
    --latency-wait 60 \
    --configfile "$CFG" \
    --printshellcmds -p
} |& tee "$LOG"

# ---- contig→bin map
BINS="$OUTDIR/bins"
MAP_RAW="$OUTDIR/contigs_to_bin.tsv"
MAP_SET="$OUTDIR/contigs_to_bin.with_set.tsv"
: > "$MAP_RAW"; shopt -s nullglob
for fa in "$BINS"/*.fa "$BINS"/*.fasta "$BINS"/*.fna; do
  bin="$(basename "${fa%.*}")"
  awk -v b="$bin" 'BEGIN{RS=">";FS="\n"} NR>1{split($1,h," "); print b "\t" h[1]}' "$fa" >> "$MAP_RAW"
done
awk '{print $1"\t"$2"\tbinny"}' "$MAP_RAW" > "$MAP_SET"

echo "[binny] $SAMPLE done → $OUTDIR"
echo "[binny] log: $LOG"
