#!/usr/bin/env bash
# geNomad filter options from config (shared by the plasmid, virus and typer steps).
#   genomad_opts  -> prints the extra end-to-end arguments, one per line
# Keys: modules.genomad.preset (default | conservative | relaxed),
#       modules.genomad.min_score (geNomad default 0.7; ignored with a preset).
# shellcheck source=bin/lib/cfg.sh
source "$ROOT/bin/lib/cfg.sh"
genomad_opts(){
  local preset score
  preset="$(cfg modules.genomad.preset default)"
  score="$(cfg modules.genomad.min_score 0.7)"
  case "$preset" in
    conservative) echo "--conservative" ;;
    relaxed)      echo "--relaxed" ;;
    default)      [[ "$score" != "0.7" ]] && printf -- '--min-score\n%s\n' "$score" ;;
    *) echo "[genomad] invalid modules.genomad.preset '$preset'" >&2; return 2 ;;
  esac
  return 0
}
