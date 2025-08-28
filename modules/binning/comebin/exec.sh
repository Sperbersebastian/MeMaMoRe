#!/usr/bin/env bash
# env: ROOT SAMPLE [THREADS] [PARAMS_YAML] [FORCE]
set -euo pipefail

ts(){ date -u +"%Y-%m-%dT%H:%M:%SZ"; }
j(){ echo "{\"ts\":\"$(ts)\",\"module\":\"binning\",\"program\":\"COMEBin\",\"sample\":\"$SAMPLE\",\"msg\":\"$1\"}"; }

# ── Inputs
ROOT="${ROOT:-/media/Box/MeMaMoRe}"
SAMPLE="${SAMPLE:?SAMPLE required}"
THREADS="${THREADS:-}"
PARAMS_YAML="${PARAMS_YAML:-}"
FORCE="${FORCE:-0}"

# ── Micromamba (pin to ~/micromamba)
MM="$HOME/micromamba/bin/micromamba"
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
command -v "$MM" >/dev/null || { echo "[COMEBin] micromamba not at $MM"; exit 3; }

# ── Threads
if [[ -z "$THREADS" ]]; then
  CFG=""
  [[ -n "$PARAMS_YAML" && -s "$PARAMS_YAML" ]] && CFG="$PARAMS_YAML"
  [[ -z "$CFG" && -s "$ROOT/config/defaults.yaml" ]] && CFG="$ROOT/config/defaults.yaml"
  [[ -z "$CFG" && -s "$ROOT/config/config.yaml"   ]] && CFG="$ROOT/config/config.yaml"
  if [[ -n "$CFG" ]]; then
    set +e
    THREADS="$("$MM" run -n env_sra_tools python "$ROOT/scripts/merge_params.py" --defaults "$CFG" --get global.threads 2>/dev/null)"
    set -e
  fi
  THREADS="${THREADS:-16}"
fi

# ── Paths
ASM_DIR="$ROOT/SRA/assemblies/spades/$SAMPLE"
CONTIGS="$ASM_DIR/contigs.fasta"
BAM_DIR="$ASM_DIR/map"
OUT_DIR="$ROOT/SRA/binning/comebin/$SAMPLE"
LOG="$OUT_DIR/run.log"
mkdir -p "$OUT_DIR"

[[ -s "$CONTIGS" ]] || { echo "[COMEBin] missing $CONTIGS"; exit 2; }
[[ -d "$BAM_DIR" ]] || { echo "[COMEBin] missing $BAM_DIR"; exit 2; }
compgen -G "$BAM_DIR"/*.bam >/dev/null || { echo "[COMEBin] no BAMs in $BAM_DIR"; exit 2; }

# ── Resolve runner: env → repo fallback
COMEBIN_DIR="$ROOT/external/comebin"
RUNNER_CMD=()
if "$MM" run -n comebin_env bash -lc 'command -v run_comebin.sh >/dev/null'; then
  RUNNER_CMD=("$MM" run -n comebin_env run_comebin.sh)
elif "$MM" run -n comebin_env bash -lc 'command -v comebin >/dev/null'; then
  RUNNER_CMD=("$MM" run -n comebin_env comebin)
elif [[ -x "$COMEBIN_DIR/run_comebin.sh" ]]; then
  RUNNER_CMD=("$COMEBIN_DIR/run_comebin.sh")
else
  echo "[COMEBin] runner not found in env or $COMEBIN_DIR" >&2
  exit 3
fi

# ── Detect existing outputs
have_bins(){
  compgen -G "$OUT_DIR/comebin_res/comebin_res_bins/"*.fa >/dev/null || \
  compgen -G "$OUT_DIR/bins/"*.fa >/dev/null
}

# ── Run if needed
if [[ "$FORCE" == "1" ]] || ! have_bins; then
  j "start"
  set -o pipefail
  "${RUNNER_CMD[@]}" \
    -a "$CONTIGS" \
    -p "$BAM_DIR" \
    -o "$OUT_DIR" \
    -t "$THREADS" | tee "$LOG"
  set +o pipefail
else
  j "bins already exist → skipping COMEBin run"
fi

# ── Export mapping for MAGScoT: contig \t bin \t set
MAP_TSV="$OUT_DIR/contigs_to_bin.with_set.tsv"
: > "$MAP_TSV"

FOUND=0
for PAT in "$OUT_DIR/comebin_res/comebin_res_bins/"*.fa "$OUT_DIR/bins/"*.fa; do
  if compgen -G "$PAT" >/dev/null; then
    FOUND=1
    for f in $PAT; do
      [[ -s "$f" ]] || continue
      bin="$(basename "${f%.*}")"
      awk '/^>/{sub(/^>/,""); sub(/ .*/,""); print}' "$f" \
        | awk -v b="$bin" 'NF{print $0 "\t" b "\tcomebin"}' >> "$MAP_TSV"
    done
  fi
done
[[ $FOUND -eq 1 ]] || echo "[COMEBin] no bin FASTAs found under comebin_res/... or bins/"

j "done → $OUT_DIR"
echo "$OUT_DIR"
