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
OUTDIR="$ROOT/SRA/binning/binny/$SAMPLE"; LOGS="$ROOT/logs"
mkdir -p "$OUTDIR" "$LOGS"

# ---- persistent rule-env prefix
PERSIST="$(abspath "$ROOT/.snk_envs/_snk_binny_persist")"
mkdir -p "$PERSIST"
export SNAKEMAKE_CONDA_PREFIX="$PERSIST"

# ---- params
CPUS="$(micromamba run -n env_binning yq -r '.modules.binning.cpus // 8' "$PARAMS_YAML" 2>/dev/null || echo 8)"
MINLEN="$(micromamba run -n env_binning yq -r '.modules.binning.min_contig_len // 1500' "$PARAMS_YAML" 2>/dev/null || echo 1500)"

# ---- prefer precomputed depth (strip header)
CDEPTH="$ROOT/SRA/binning/coverage/$SAMPLE/$SAMPLE.depth.txt"
CDEPTH_NOHDR="$ROOT/SRA/binning/coverage/$SAMPLE/$SAMPLE.depth.nohdr.txt"
USE_DEPTH=0
if [[ -s "$CDEPTH" ]]; then
  if [[ ! -s "$CDEPTH_NOHDR" ]]; then
    tail -n +2 "$CDEPTH" > "$CDEPTH_NOHDR"
  fi
  USE_DEPTH=1
fi

# ---- MANTIS custom-only wiring (local assets only)
export NLTK_DATA="$BINNY_DIR/database/support/nltk_data"
cfg_mantis="$BINNY_DIR/config/binny_mantis.cfg"
sed -i -E \
  -e "s|^custom_ref=.*/checkm_tf/checkm_filtered_tf\.hmm$|custom_ref=$BINNY_DIR/database/hmms/checkm_tf/checkm_filtered_tf.hmm|" \
  -e "s|^custom_ref=.*/checkm_pf/checkm_filtered_pf\.hmm$|custom_ref=$BINNY_DIR/database/hmms/checkm_pf/checkm_filtered_pf.hmm|" \
  "$cfg_mantis" || true
grep -q '^checkm_filtered_tf_weight=' "$cfg_mantis" || echo "checkm_filtered_tf_weight=0.5" >>"$cfg_mantis"
grep -q '^checkm_filtered_pf_weight=' "$cfg_mantis" || echo "checkm_filtered_pf_weight=1"   >>"$cfg_mantis"
echo "[binny] MANTIS custom-only mode (no DB downloads)."

# ---- Binny config.yaml for Snakefile
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

# ---- make helpers importable and provide conda CLI
export PYTHONPATH="$BINNY_DIR:$BINNY_DIR/workflow/scripts:${PYTHONPATH:-}"
export CONDA_EXE="$(micromamba run -n env_binning bash -lc 'command -v conda')"
mkdir -p "$BINNY_DIR/logs"


# ---- run Snakemake and tee logs (stdout+stderr)
{
  snakemake -s "$BINNY_DIR/Snakefile" \
    --directory "$BINNY_DIR" \
    --cores "$CPUS" \
    --use-conda \
    --conda-prefix "$SNAKEMAKE_CONDA_PREFIX" \
    --configfile "$CFG" \
    --printshellcmds -p
} |& tee "$LOGS/binning_binny_${SAMPLE}.log"

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
echo "[binny] log: $LOGS/binning_binny_${SAMPLE}.log"
