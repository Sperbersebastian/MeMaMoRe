#!/usr/bin/env bash
# cluster/exec.sh — global clustering of all per-sample plasmid representatives
set -euo pipefail
ROOT="${1:?}"; SAMPLE_IGNORED="${2:-_}"; CPUS="${3:-16}"; FORCE="${4:-0}"

GOUT="$ROOT/SRA/plasmids/_global/cluster"; mkdir -p "$GOUT"
GCONCAT="$GOUT/plasmids_concat.all_samples.fasta"
GREP_OUT="$GOUT/plasmids_derep.global.fasta"

ID="${ID:-95}"
COV="${COV:-80}"

CLSTR="${CLSTR:-$(find "$ROOT/external/Stampede-ClusterGenomes" -type f \( -name 'Cluster_genomes.pl' -o -path '*Cluster_genomes.pl' \) -print -quit)}"
[[ -n "${CLSTR:-}" && -x "$CLSTR" ]] || { echo "[err] Cluster_genomes.pl not found/executable"; exit 2; }

# 1) All-sample set, rebuilt from every sample's dereplicated plasmids. Headers
#    get the sample as prefix ("<sample>|<tool>|<contig>"): contig names repeat
#    across samples, and the sample is needed to link back (HOTSPOT, hosts).
tmp="$GOUT/.tmp.concat.$$.fa"; : > "$tmp"
shopt -s nullglob
for rep in "$ROOT"/SRA/plasmids/*/cluster/plasmids_derep.*.fasta; do
  [[ -s "$rep" ]] || continue
  s="$(basename "$(dirname "$(dirname "$rep")")")"
  [[ "$s" == _* ]] && continue
  awk -v s="$s" '/^>/{sub(/^>/, ">" s "|")} {print}' "$rep" >> "$tmp"
done
shopt -u nullglob
[[ -s "$tmp" ]] || { rm -f "$tmp"; echo "[global] no per-sample plasmid representatives yet"; exit 0; }
if [[ -s "$GCONCAT" ]] && cmp -s "$tmp" "$GCONCAT"; then
  rm -f "$tmp"
else
  mv -f "$tmp" "$GCONCAT"
  # set changed: recluster. Stampede reuses an existing <stem>-nucmer.out.* /
  # -cover.csv without checking the input, so those must go as well.
  rm -f "$GOUT"/.done.global.* "$GOUT"/plasmids_concat*_"${ID}"-"${COV}".{fna,clstr} \
        "$GOUT"/plasmids_concat.all_samples-nucmer.out.* "$GOUT"/plasmids_concat.all_samples-cover.csv
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
# Stampede writes <input stem>_<ID>-<COV>.fna and .clstr (input plasmids_concat.all_samples.fasta)
GCENTROID="$(ls -1 "$GOUT"/plasmids_concat*_"${ID}"-"${COV}".fna 2>/dev/null | head -n1 || true)"
if [[ -n "${GCENTROID:-}" && -s "$GCENTROID" ]]; then
  cp -f "$GCENTROID" "$GREP_OUT"
else
  GREP_IDS="$GOUT/global.rep.ids"
  # representatives = 2nd field of the ">Cluster_N" lines of the .clstr file
  # (global.clusters.txt is Stampede's verbose log, not a cluster list)
  GCLSTR="$(ls -1 "$GOUT"/plasmids_concat*_"${ID}"-"${COV}".clstr 2>/dev/null | head -n1 || true)"
  [[ -n "$GCLSTR" ]] || { echo "[err] no Stampede .clstr in $GOUT"; exit 4; }
  awk -F'\t' '/^>Cluster/ && NF>=2 {print $2}' "$GCLSTR" > "$GREP_IDS"

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
