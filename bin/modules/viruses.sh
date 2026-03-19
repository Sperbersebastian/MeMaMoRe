#!/usr/bin/env bash
# Virus Detection Module Coordinator (bin/modules/viruses.sh)
# Per-sample pipeline: metaviralSPAdes → viralVerify → geNomad → VirSorter2 → VIBRANT → union → CheckV
set -euo pipefail

: "${ROOT:?ROOT not set}"
: "${SAMPLE:?SAMPLE not set}"
CPUS="${CPUS:-16}"
FORCE="${FORCE:-0}"

# ---------- micromamba ----------
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
command -v micromamba >/dev/null 2>&1 || { echo "[error] micromamba not found"; exit 1; }
eval "$(micromamba shell hook --shell=bash)"

# env manager
# shellcheck source=bin/lib/env.sh
source "$ROOT/bin/lib/env.sh"

# ---------- helpers ----------
_run_step(){
  local step="$1" exec_path="$2" sample="$3" ; shift 3
  local OUT="$ROOT/SRA/viruses/$sample"; mkdir -p "$OUT/.done"
  local done="$OUT/.done/${step}.done"

  if [[ -s "$done" && "$FORCE" != "1" ]]; then
    echo "[skip] $sample :: $step (done)"; return 0
  fi
  ensure_env_by_module "$step"
  echo "[run]  $sample :: $step -> $exec_path"
  local t0; t0=$(date +%s)
  local env_name
  case "$step" in
    metaviralspades)  env_name="env_assembly_core" ;;
    viralverify)      env_name="viralverify_env" ;;
    genomad)          env_name="genomad_env" ;;
    virsorter2)       env_name="virsorter2_env" ;;
    vibrant)          env_name="vibrant_env" ;;
    union_virus)      env_name="plasmids_core" ;;
    checkv)           env_name="checkv_env" ;;
    *)                env_name="plasmids_core" ;;
  esac
  run_in_env "$env_name" /bin/bash "$exec_path" "$ROOT" "$sample" "$CPUS" "$FORCE"
  date -u +"%Y-%m-%dT%H:%M:%SZ" >"$done"
  echo "[ok]   $sample :: $step ($(( $(date +%s) - t0 ))s)"
}

# ---------- logging ----------
MOD="$ROOT/modules/viruses"
ASM_MOD="$ROOT/modules/viruses/metaviralspades/exec.sh"
OUTDIR="$ROOT/SRA/viruses/$SAMPLE"; mkdir -p "$OUTDIR/.done"
LOGDIR="$ROOT/logs"; mkdir -p "$LOGDIR"
LOGFILE="$LOGDIR/viruses_${SAMPLE}.log"
exec > >(stdbuf -oL tee -a "$LOGFILE") 2>&1
echo "[viruses] START sample=$SAMPLE cpus=$CPUS force=$FORCE host=$(hostname)"

# verify scripts exist
[[ -x "$ASM_MOD" ]] || { echo "[error] missing: $ASM_MOD"; exit 1; }
for s in viralverify genomad virsorter2 vibrant union checkv; do
  [[ -x "$MOD/$s/exec.sh" ]] || { echo "[error] missing: $MOD/$s/exec.sh"; exit 1; }
done

# ---------- execution order ----------
# 0) metaviralSPAdes
_run_step "metaviralspades" "$ASM_MOD" "$SAMPLE"
# 1) viralVerify on metaviralSPAdes contigs
_run_step "viralverify"     "$MOD/viralverify/exec.sh"  "$SAMPLE"
# 2) geNomad
_run_step "genomad"         "$MOD/genomad/exec.sh"      "$SAMPLE"
# 3) VirSorter2
_run_step "virsorter2"      "$MOD/virsorter2/exec.sh"   "$SAMPLE"
# 4) VIBRANT
_run_step "vibrant"         "$MOD/vibrant/exec.sh"      "$SAMPLE"
# 5) union + dedup
_run_step "union_virus"     "$MOD/union/exec.sh"        "$SAMPLE"
# 6) CheckV quality
_run_step "checkv"          "$MOD/checkv/exec.sh"       "$SAMPLE"

echo "[viruses] DONE sample=$SAMPLE"
