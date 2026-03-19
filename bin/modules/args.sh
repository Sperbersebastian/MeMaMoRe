#!/usr/bin/env bash
# ARG Detection Module Coordinator (bin/modules/args.sh)
# Per-sample pipeline: DeepARG → RGI → AMR++ → summary
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
  local OUT="$ROOT/SRA/args/$sample"; mkdir -p "$OUT/.done"
  local done="$OUT/.done/${step}.done"

  if [[ -s "$done" && "$FORCE" != "1" ]]; then
    echo "[skip] $sample :: $step (done)"; return 0
  fi
  ensure_env_by_module "$step"
  echo "[run]  $sample :: $step -> $exec_path"
  local t0; t0=$(date +%s)
  local env_name
  case "$step" in
    deeparg)       env_name="deeparg_env" ;;
    rgi)           env_name="rgi_env" ;;
    amrplusplus)   env_name="amrplusplus_env" ;;
    summary_args)  env_name="deeparg_env" ;;
    *)             env_name="deeparg_env" ;;
  esac
  run_in_env "$env_name" /bin/bash "$exec_path" "$ROOT" "$sample" "$CPUS" "$FORCE"
  date -u +"%Y-%m-%dT%H:%M:%SZ" >"$done"
  echo "[ok]   $sample :: $step ($(( $(date +%s) - t0 ))s)"
}

# ---------- logging ----------
MOD="$ROOT/modules/args"
OUTDIR="$ROOT/SRA/args/$SAMPLE"; mkdir -p "$OUTDIR/.done"
LOGDIR="$ROOT/logs"; mkdir -p "$LOGDIR"
LOGFILE="$LOGDIR/args_${SAMPLE}.log"
exec > >(stdbuf -oL tee -a "$LOGFILE") 2>&1
echo "[args] START sample=$SAMPLE cpus=$CPUS force=$FORCE host=$(hostname)"

# verify scripts exist
for s in deeparg rgi amrplusplus summary; do
  [[ -x "$MOD/$s/exec.sh" ]] || { echo "[error] missing: $MOD/$s/exec.sh"; exit 1; }
done

# ---------- execution order ----------
# 0) DeepARG
_run_step "deeparg"       "$MOD/deeparg/exec.sh"       "$SAMPLE"
# 1) RGI (CARD)
_run_step "rgi"           "$MOD/rgi/exec.sh"           "$SAMPLE"
# 2) AMR++
_run_step "amrplusplus"   "$MOD/amrplusplus/exec.sh"   "$SAMPLE"
# 3) Summary
_run_step "summary_args"  "$MOD/summary/exec.sh"       "$SAMPLE"

echo "[args] DONE sample=$SAMPLE"
