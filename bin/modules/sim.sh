#!/usr/bin/env bash
set -euo pipefail
: "${ROOT:?ROOT not set}"
CFG="${CFG:-$ROOT/config/sim.yaml}"

# use env helpers
# shellcheck source=bin/lib/env.sh
source "$ROOT/bin/lib/env.sh"

# ensure sim env exists
ensure_env_by_module "sim"

# wrapper to run yq inside sim_env
yq_in(){ run_in_env "sim_env" yq -r "$1" "$CFG"; }

# parse config via yq in env
SET_NAME="$(yq_in '.set_name')"
N_READS="$(yq_in '.n_reads')"
READ_LEN="$(yq_in '.read_len')"
MODEL="$(yq_in '.model')"
CPUS="$(yq_in '.cpus')"
SEED="$(yq_in '.seed')"

HOST_HQ="$(yq_in '.hosts.hq')"
HOST_LQ="$(yq_in '.hosts.lq')"
PLASMIDS="$(yq_in '.plasmids[]' | xargs)"
PHAGES="$(yq_in '.phages[]' | xargs)"

# export abundances as ABUND_<accession-with-dots-replaced>
while IFS=$'\t' read -r acc w; do
  [[ -z "$acc" || "$acc" == "null" ]] && continue
  safe=${acc//./_}
  export "ABUND_${safe}=$w"
done < <(run_in_env "sim_env" yq -r '.abundances | to_entries[] | "\(.key)\t\(.value)"' "$CFG")

export SET_NAME N_READS READ_LEN MODEL CPUS SEED HOST_HQ HOST_LQ PLASMIDS PHAGES
modules/sim/exec.sh
