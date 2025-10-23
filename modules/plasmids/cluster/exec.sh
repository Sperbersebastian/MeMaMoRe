#!/usr/bin/env bash
# cluster/exec.sh — global clustering of all per-sample reps; builds concat if missing
set -euo pipefail
ROOT="${1:?}"; SAMPLE_IGNORED="${2:-_}"; CPUS="${3:-16}"; FORCE="${4:-0}"

GOUT="$ROOT/SRA/plasmids/_global/cluster"; mkdir -p "$GOUT"
GCONCAT="$GOUT/plasmids_concat.all_samples.fasta"
GREP_OUT="$GOUT/plasmids_derep.global.fasta"

ID="${ID:-95}"
COV="${COV:-80}"

CLSTR="${CLSTR:-$(find "$ROOT/external/Stampede-ClusterGenomes" -type f \( -name 'Cluster_genomes.pl' -o -path '*Cluster_genomes.pl' \) -print -quit)}"
[[ -n "${CLSTR:-}" && -x "$CLSTR" ]] || { echo "[err] Cluster_genomes.pl not found/executable"; exit 2; }

# 1) Build concat if missing by scanning per-sample reps
if [[ ! -s "$GCONCAT" ]]; then
  tmp="$GOUT/.tmp.concat.$$.fa"; : > "$tmp"
  shopt -s nullglob
  for rep in "$ROOT"/SRA/plasmids/*/cluster/plasmids_derep.*.fasta; do
    [[ -s "$rep" ]] || continue
    awk '
      /^>/ {h=$0; if(!(h in seen)){seen[h]=1; keep=1; print; next} keep=0}
      { if(keep) print }
    ' "$rep" >> "$tmp"
  done
  shopt -u nullglob
  [[ -s "$tmp" ]] || { echo "[global] keine Sample-Reps gefunden"; exit 3; }
  mv -f "$tmp" "$GCONCAT"
fi

NALL=$(grep -c '^>' "$GCONCAT" || true)
echo "[global] concat ready: $NALL seq -> $GCONCAT"

# 2) Global Stampede clustering
GFLAG="$GOUT/.done.global.i${ID}.c${COV}"
[[ "$FORCE" == "1" ]] && rm -f "$GFLAG" "$GOUT"/plasmids_concat_*"${ID}"-"${COV}".fna

if [[ ! -f "$GFLAG" ]]; then
  echo "[global] clustering id=$ID cov=$COV"
  ( cd "$GOUT" && perl "$CLSTR" -f "plasmids_concat.all_samples.fasta" -i "$ID" -c "$COV" > global.clusters.txt 2> global.stderr )
  touch "$GFLAG"
else
  echo "[global] cached -> $GOUT/global.clusters.txt"
fi
[[ -s "$GOUT/global.clusters.txt" ]] || { echo "[err] global.clusters.txt leer"; exit 4; }

# 3) Representatives (prefer Stampede centroid)
GCENTROID="$(ls -1 "$GOUT"/plasmids_concat_*"${ID}"-"${COV}".fna 2>/dev/null | head -n1 || true)"
if [[ -n "${GCENTROID:-}" && -s "$GCENTROID" ]]; then
  cp -f "$GCENTROID" "$GREP_OUT"
else
  GREP_IDS="$GOUT/global.rep.ids"
  awk '
    BEGIN{picked=0}
    /^[[:space:]]*$/ {next}
    /^#/ {next}
    /^[>]*[Cc]luster[[:space:]]*[0-9]+/ {picked=0; next}
    { if(!picked){ id=$1; sub(/^[>]/,"",id); print id; picked=1 } }
  ' "$GOUT/global.clusters.txt" > "$GREP_IDS"

  awk 'BEGIN{
         while((getline k<ARGV[1])>0){want[k]=1}
         close(ARGV[1]); ARGV[1]=""
       }
       /^>/{
         id=substr($0,2); sub(/[ \t].*$/,"",id)
         keep=(id in want)
       }
       { if(keep) print }
  ' "$GREP_IDS" "$GCONCAT" > "$GREP_OUT"
fi

NUNIQ=$(grep -c '^>' "$GREP_OUT" || true)
echo "[global-derep] $NUNIQ / $NALL reps -> $GREP_OUT"
