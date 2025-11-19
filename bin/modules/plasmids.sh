#!/usr/bin/env bash
# Plasmid-Workflow-Orchestrator (bin/modules/plasmids.sh)
# Modes:
#   per-sample (default): requires ROOT, SAMPLE
#   --global-cluster    : requires ROOT; builds global derep from concat and optional per-sample CoverM
set -euo pipefail

# ---------- arg parsing ----------
MODE="sample"
if [[ "${1:-}" == "--global-cluster" ]]; then
  MODE="global"; shift
fi

: "${ROOT:?ROOT not set}"
if [[ "$MODE" == "sample" ]]; then
  : "${SAMPLE:?SAMPLE not set}"
fi
CPUS="${CPUS:-16}"
FORCE="${FORCE:-0}"
FORCE_ENV="${FORCE_ENV:-${force_env:-${forece_env:-}}}"
RUN_COVERM_AFTER_GLOBAL="${RUN_COVERM_AFTER_GLOBAL:-0}"
SAMPLES_FILE="${SAMPLES_FILE:-}"

# ---------- micromamba ----------
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
command -v micromamba >/dev/null 2>&1 || { echo "[error] micromamba not found"; exit 1; }
eval "$(micromamba shell hook --shell=bash)"

if [[ -n "${FORCE_ENV:-}" ]]; then
  export FORCE_ENV_MODE="one"
  export FORCE_ENV_TARGET="$FORCE_ENV"
fi

# env manager
# shellcheck source=bin/lib/env.sh
source "$ROOT/bin/lib/env.sh"

# ---------- helpers ----------
_have_reads(){
  compgen -G "$ROOT/SRA/qc/fastp/$1/*_R1*.fastq.gz" >/dev/null || \
  compgen -G "$ROOT/SRA/qc/fastp/$1/*trimmed_1*.fastq.gz" >/dev/null
}



_need_bin(){
  case "$1" in
    metaplasmidspades) echo ""metaplasmidspades.py"" ;;
    viralverify|plasme|typing_hotspot) echo "python" ;;
    union_cluster|cluster|global_cluster) echo "perl" ;;
    genomad)     echo "genomad" ;;
    mobrecon)    echo "mob_recon" ;;
    cluster)     echo "perl" ;;
    coverm)      echo "coverm" ;;
    typer)       echo "mob_typer" ;;
    *)           echo "" ;;
  esac
}

_env_for(){
  case "$1" in
    metaplasmidspades)      echo "env_assembly_core" ;;
    viralverify)            echo "viralverify_env" ;;
    plasme)                 echo "plasme_env" ;;
    genomad)                echo "genomad_env" ;;
    mobrecon|typer)         echo "mobsuite_env" ;;
    cluster|global_cluster|union_cluster) echo "stampede_env" ;;
    coverm|union_cluster)   echo "plasmids_core" ;;
    typing_hotspot)         echo "hotspot_env" ;;
    *)                      echo "plasmids_core" ;;
  esac
}

# Cache for environment and binary checks to avoid repeated lookups
declare -A _ENV_CACHE _BIN_CACHE

_run_step(){
  local step="$1" exec_path="$2" sample="$3" ; shift 3
  local OUT="$ROOT/SRA/plasmids/$sample"; mkdir -p "$OUT/.done"
  local done="$OUT/.done/${step}.done"

  if [[ -s "$done" && "$FORCE" != "1" ]]; then
    echo "[skip] $sample :: $step (done)"; return 0
  fi
  
  # Use cached environment name or compute and cache it
  local env_name="${_ENV_CACHE[$step]:-}"
  if [[ -z "$env_name" ]]; then
    ensure_env_by_module "$step"
    env_name="$(_env_for "$step")"
    _ENV_CACHE[$step]="$env_name"
  fi
  
  # Use cached binary check or compute and cache it
  local need="${_BIN_CACHE[$step]:-NOT_SET}"
  if [[ "$need" == "NOT_SET" ]]; then
    need="$(_need_bin "$step")"
    _BIN_CACHE[$step]="$need"
  fi
  
  if [[ -n "$need" ]]; then
    run_in_env "$env_name" bash -lc "command -v $need >/dev/null" \
      || { echo "[error] '$need' missing in env '$env_name' (step $step)"; exit 127; }
  fi
  echo "[run]  $sample :: $step -> $exec_path $*"
  local t0; t0=$(date +%s)
  run_in_env "$env_name" /bin/bash "$exec_path" "$ROOT" "$sample" "$CPUS" "$FORCE"
  date -u +"%Y-%m-%dT%H:%M:%SZ" >"$done"
  echo "[ok]   $sample :: $step ($(( $(date +%s) - t0 ))s)"
}

_global_cluster(){
  ensure_env_by_module cluster
  local GDIR="$ROOT/SRA/plasmids/_global/cluster"
  local IN="$GDIR/plasmids_concat.all_samples.fasta"
  local OUT="$GDIR"
  mkdir -p "$GDIR"
  [[ -s "$IN" ]] || { echo "[global_cluster] no concat: $IN"; exit 1; }

  local ID="${ID:-95}"
  local COV="${COV:-80}"
  local ST_CL="$ROOT/external/Stampede-ClusterGenomes/Cluster_genomes.pl"
  [[ -x "$ST_CL" ]] || { echo "[global_cluster] Cluster_genomes.pl missing at $ST_CL"; exit 2; }

  run_in_env stampede_env bash -lc \
    "perl '$ST_CL' -f '$IN' -i '$ID' -c '$COV' >'$OUT/global.clusters.txt' 2>'$OUT/global.stderr'"

  [[ -s "$OUT/global.clusters.txt" ]] || { echo "[global_cluster] empty clusters"; exit 3; }

  # Optimized: combine extraction and filtering into a single AWK pass
  awk -v infile="$IN" '
    # First pass: extract rep IDs from clusters.txt
    FILENAME==ARGV[1] {
      if (/^[[:space:]]*$/ || /^#/) next
      if (/^[>]*[[:space:]]*[Cc]luster[[:space:]]*[0-9]+/) {picked=0; next}
      if (!picked) {
        id=$1; sub(/^[>]/,"",id)
        want[id]=1
        picked=1
      }
      next
    }
    # Second pass: filter sequences from input FASTA
    /^>/ {
      id=substr($0,2); sub(/[ \t].*$/,"",id)
      keep=(id in want)
    }
    { if(keep) print }
  ' "$OUT/global.clusters.txt" "$IN" > "$OUT/plasmids_derep.global.fasta"

  local NALL; NALL=$(grep -c '^>' "$IN" || true)
  local NUNQ; NUNQ=$(grep -c '^>' "$OUT/plasmids_derep.global.fasta" || true)
  echo "[global_cluster] unique $NUNQ / $NALL -> $OUT/plasmids_derep.global.fasta"

  if [[ "$RUN_COVERM_AFTER_GLOBAL" == "1" ]]; then
    ensure_env_by_module coverm
    local FA="$OUT/plasmids_derep.global.fasta"
    [[ -s "$FA" ]] || { echo "[coverm] missing global fasta: $FA"; exit 2; }

    local samples=()
    if [[ -n "$SAMPLES_FILE" && -r "$SAMPLES_FILE" ]]; then
      mapfile -t samples < <(grep -v '^\s*$' "$SAMPLES_FILE")
    else
      samples=( $(find "$ROOT/SRA/plasmids" -maxdepth 2 -type d -name cluster -printf '%h\n' \
                  | awk -F/ '{print $(NF)}' | grep -v '^_global$' | sort -u) )
    fi

    for s in "${samples[@]}"; do
      if _have_reads "$s"; then
        echo "[coverm] $s"
        [[ -x "$ROOT/modules/plasmids/coverm/exec.sh" ]] || { echo "[error] missing: modules/plasmids/coverm/exec.sh"; exit 1; }
        _run_step "coverm" "$ROOT/modules/plasmids/coverm/exec.sh" "$s"
      else
        echo "[coverm] skip $s (no reads)"
      fi
    done
  fi
}

# ---------- execution ----------
if [[ "$MODE" == "global" ]]; then
  echo "[plasmids] GLOBAL-CLUSTER cpus=$CPUS"
  _global_cluster
  exit 0
fi

# ensure assembly core (spades + metaPlasmidSPAdes)
ensure_env_by_module assembly
run_in_env env_assembly_core bash -lc 'command -v "metaplasmidspades.py" >/dev/null || { echo "[error] metaplasmidspades.py missing in env_assembly_core"; exit 2; }'

# logging
MOD="$ROOT/modules/plasmids"
ASM_MOD="$ROOT/modules/assembly/metaplasmidspades/exec.sh"
OUTDIR="$ROOT/SRA/plasmids/$SAMPLE"; mkdir -p "$OUTDIR/.done"
LOGDIR="$ROOT/logs"; mkdir -p "$LOGDIR"
LOGFILE="$LOGDIR/plasmids_${SAMPLE}.log"
exec > >(stdbuf -oL tee -a "$LOGFILE") 2>&1
echo "[plasmids] START sample=$SAMPLE cpus=$CPUS force=$FORCE host=$(hostname)"

# verify scripts
[[ -x "$ASM_MOD" ]] || { echo "[error] missing: modules/assembly/metaplasmidspades/exec.sh"; exit 1; }
for s in viralverify plasme genomad mobrecon union_cluster cluster typing_hotspot typer; do
  [[ -x "$MOD/$s/exec.sh" ]] || { echo "[error] missing: $MOD/$s/exec.sh"; exit 1; }
done
[[ -x "$MOD/coverm/exec.sh" ]] || true

# ORDER per sample:
# 0) metaPlasmidSPAdes assembly for VV input
_run_step "metaplasmidspades" "$ASM_MOD" "$SAMPLE"
# 1) viralVerify on metaPlasmidSPAdes contigs
_run_step "viralverify"    "$MOD/viralverify/exec.sh"    "$SAMPLE"
# 2) geNomad, PLASMe, mob_recon on metaSPAdes contigs
_run_step "genomad"        "$MOD/genomad/exec.sh"        "$SAMPLE"
_run_step "plasme"         "$MOD/plasme/exec.sh"         "$SAMPLE"
_run_step "mobrecon"       "$MOD/mobrecon/exec.sh"       "$SAMPLE"
# 3) collect plasmids from 0–2
_run_step "union_cluster"  "$MOD/union_cluster/exec.sh"  "$SAMPLE"
# 4) per-sample typing
_run_step "typer"          "$MOD/typer/exec.sh"          "$SAMPLE"
_run_step "typing_hotspot" "$MOD/typing_hotspot/exec.sh" "$SAMPLE"

echo "[plasmids] DONE sample=$SAMPLE (global cluster: bin/modules/plasmids.sh --global-cluster ROOT=/path)"
