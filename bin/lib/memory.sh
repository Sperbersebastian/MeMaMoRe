#!/usr/bin/env bash
# Memory gate for RAM-hungry steps (GTDB-Tk, Plasmer).
#   mem_gate <tool> <min_gb>   -> 0 = run, 1 = skip (message printed)
# Counts RAM + swap: on a 15 GB laptop with 50 GB swap both tools finished
# (peak swap 6-9 GB), only slower. <TOOL>_FORCE=1 runs regardless.
mem_gate(){
  local tool="$1" min="$2" ram swap force
  ram=$(( $(awk '/^MemTotal:/{print $2}' /proc/meminfo 2>/dev/null || echo 0) / 1024 / 1024 ))
  swap=$(( $(awk '/^SwapTotal:/{print $2}' /proc/meminfo 2>/dev/null || echo 0) / 1024 / 1024 ))
  force="$(echo "${tool}_FORCE" | tr '[:lower:]' '[:upper:]')"
  if [[ "${!force:-0}" == "1" ]]; then return 0; fi
  if (( ram + swap < min )); then
    echo "[$tool] SKIP: ${ram} GB RAM + ${swap} GB swap < ${min} GB needed (${force}=1 overrides)"
    return 1
  fi
  if (( ram < min )); then
    echo "[$tool] note: ${ram} GB RAM < ${min} GB, relying on ${swap} GB swap (slower)"
  fi
  return 0
}
