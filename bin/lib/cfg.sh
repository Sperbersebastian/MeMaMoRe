#!/usr/bin/env bash
# Shared config lookup for module scripts.
#   cfg <dotted.key> <default>
# Reads the merged PARAMS_YAML (config/defaults.yaml + config/local.yaml + --set,
# built by bin/main.sh); falls back to <default> when the key is missing or the
# script runs standalone without PARAMS_YAML.
cfg(){
  local v=""
  if [[ -n "${PARAMS_YAML:-}" && -s "$PARAMS_YAML" ]]; then
    v="$(micromamba run -n env_sra_tools python "$ROOT/scripts/merge_params.py" \
           --defaults "$PARAMS_YAML" --get "$1" 2>/dev/null || true)"
  fi
  if [[ -n "$v" && "$v" != "None" && "$v" != "null" ]]; then echo "$v"; else echo "$2"; fi
}
