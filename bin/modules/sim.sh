#!/usr/bin/env bash
set -euo pipefail
: "${ROOT:?ROOT not set}"
CFG="${CFG:-$ROOT/config/sim.yaml}"

# use env helpers
# shellcheck source=bin/lib/env.sh
source "$ROOT/bin/lib/env.sh"

# ensure sim env exists
ensure_env_by_module "sim"

# config -> plain tables (references, per-sample weights, truth); see scripts/sim_plan.py
PLAN="$(mktemp -d)"
trap 'rm -rf "$PLAN"' EXIT
run_in_env "sim_env" python "$ROOT/scripts/sim_plan.py" "$CFG" "$PLAN"

setting(){ awk -F'\t' -v k="$1" '$1==k{print $2}' "$PLAN/settings.tsv"; }
SET_NAME="$(setting set_name)"
N_READS="$(setting n_reads)"
READ_LEN="$(setting read_len)"
SEED="$(setting seed)"
MODEL="$(setting model)"
CPUS="$(setting cpus)"

export SET_NAME N_READS READ_LEN MODEL CPUS SEED PLAN
bash "$ROOT/modules/sim/exec.sh"
