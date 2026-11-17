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

#   genomad_restart IN OUT FORCE -> prints "--restart" when FORCE=1 or when OUT holds
#   intermediate files older than IN (new contigs): geNomad otherwise skips every step
#   whose output exists ("Previous execution detected") and reports on stale data.
genomad_restart(){
  local in="$1" out="$2" force="${3:-0}"
  if [[ "$force" == "1" ]] || [[ -n "$(find "$out" -mindepth 2 -type f ! -newer "$in" -print -quit 2>/dev/null)" ]]; then
    echo "--restart"
  fi
  return 0
}
