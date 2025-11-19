#!/usr/bin/env bash
# union_cluster/exec.sh — per-sample union + clustering, then update global and trigger global clustering
set -euo pipefail
ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

BASE="$ROOT/SRA/plasmids/$SAMPLE"
UOUT="$BASE/union";   mkdir -p "$UOUT"
SOUT="$BASE/cluster"; mkdir -p "$SOUT"
GOUT="$ROOT/SRA/plasmids/_global/cluster"; mkdir -p "$GOUT"

RAW="$UOUT/plasmids_raw.fasta"
LRAW="$SOUT/plasmids_raw.fasta"
REP_OUT="$SOUT/plasmids_derep.$SAMPLE.fasta"
GCONCAT="$GOUT/plasmids_concat.all_samples.fasta"

MIN_LEN="${MIN_LEN:-2000}"
ID="${ID:-95}"
COV="${COV:-80}"

# 1) Build union (no Stampede needed)
[[ "$FORCE" == "1" ]] && rm -f "$RAW"
if [[ ! -s "$RAW" ]]; then
  tmp="$UOUT/.gather.$$.fa"
  
  # Combine all sources with their prefixes into a single temp file, then filter by length
  {
    # viralVerify
    [[ -s "$BASE/viralverify/Prediction_results_fasta/contigs_plasmid.fasta" ]] && \
      awk -v p="vv" '/^>/{sub(/^>/,">"p"|")}1' \
         "$BASE/viralverify/Prediction_results_fasta/contigs_plasmid.fasta"
    [[ -s "$BASE/viralverify/Prediction_results_fasta/contigs_plasmid_uncertain.fasta" ]] && \
      awk -v p="vv_uncertain" '/^>/{sub(/^>/,">"p"|")}1' \
         "$BASE/viralverify/Prediction_results_fasta/contigs_plasmid_uncertain.fasta"

    # geNomad
    [[ -s "$BASE/genomad/contigs_summary/contigs_plasmid.fna" ]] && \
      awk -v p="genomad" '/^>/{sub(/^>/,">"p"|")}1' \
         "$BASE/genomad/contigs_summary/contigs_plasmid.fna"

    # PLASMe (prefer filtered)
    if [[ -s "$BASE/plasme/filtered_hp/plasme.hp.filtered.fasta" ]]; then
      awk -v p="plasme" '/^>/{sub(/^>/,">"p"|")}1' \
        "$BASE/plasme/filtered_hp/plasme.hp.filtered.fasta"
    else
      pf="$(ls -1 "$BASE/plasme/"*.fa "$BASE/plasme/"*.fna "$BASE/plasme/"*.fasta 2>/dev/null | head -n1 || true)"
      [[ -n "${pf:-}" && -s "$pf" ]] && awk -v p="plasme" '/^>/{sub(/^>/,">"p"|")}1' "$pf"
    fi

    # MOB-recon
    for f in "$BASE"/mobrecon/*reconstructed*.fasta "$BASE"/mobrecon/*plasmid*.fasta; do
      [[ -s "$f" ]] && awk -v p="mobrecon" '/^>/{sub(/^>/,">"p"|")}1' "$f"
    done
  } | awk -v MIN="$MIN_LEN" '
    /^>/{
      if (seqlen>=MIN && hdr!=""){print hdr; print seq}
      hdr=$0; seq=""; seqlen=0; next
    }
    {seqlen+=length($0); seq=seq $0}
    END{ if (seqlen>=MIN && hdr!=""){print hdr; print seq} }
  ' > "$RAW"

  [[ -s "$RAW" ]] || { echo "[union] no sequences after filter"; exit 3; }
  echo "[union] $(grep -c '^>' "$RAW") seq -> $RAW"
else
  echo "[union] cached -> $RAW"
fi

# 2) Per-sample clustering (Stampede)
CLSTR="${CLSTR:-$(find "$ROOT/external/Stampede-ClusterGenomes" -type f \( -name 'Cluster_genomes.pl' -o -path '*/Cluster_genomes.pl' \) -print -quit)}"
[[ -n "${CLSTR:-}" && -x "$CLSTR" ]] || { echo "[err] Cluster_genomes.pl not found/executable"; exit 2; }

cp -f "$RAW" "$LRAW"
FLAG="$SOUT/.done.sample.i${ID}.c${COV}"
[[ "$FORCE" == "1" ]] && rm -f "$FLAG" "$SOUT"/plasmids_raw_*"${ID}"-"${COV}".fna

if [[ ! -f "$FLAG" ]]; then
  echo "[stampede] $SAMPLE  id=$ID cov=$COV"
  ( cd "$SOUT" && perl "$CLSTR" -f "plasmids_raw.fasta" -i "$ID" -c "$COV" > clusters.txt 2> cluster.stderr )
  touch "$FLAG"
else
  echo "[stampede] cached -> $SOUT/clusters.txt"
fi
[[ -s "$SOUT/clusters.txt" ]] || { echo "[err] clusters.txt empty"; exit 4; }

# representatives
CENTROID_FA="$(ls -1 "$SOUT"/plasmids_raw_*"${ID}"-"${COV}".fna 2>/dev/null | head -n1 || true)"
if [[ -n "${CENTROID_FA:-}" && -s "$CENTROID_FA" ]]; then
  cp -f "$CENTROID_FA" "$REP_OUT"
else
  # Optimized: combine the two AWK passes into one by doing extraction and filtering in a single pass
  awk -v lraw="$LRAW" '
    # First pass: extract rep IDs from clusters.txt
    FILENAME==ARGV[1] {
      if (/^[[:space:]]*$/ || /^#/) next
      if (/^[>]*[Cc]luster[[:space:]]*[0-9]+/) {picked=0; next}
      if (!picked) {
        id=$1; sub(/^[>]/,"",id)
        want[id]=1
        picked=1
      }
      next
    }
    # Second pass: filter sequences from LRAW
    /^>/ {
      id=substr($0,2); sub(/[ \t].*$/,"",id)
      keep=(id in want)
    }
    { if(keep) print }
  ' "$SOUT/clusters.txt" "$LRAW" > "$REP_OUT"
fi
[[ -s "$REP_OUT" ]] || { echo "[err] no representatives"; exit 5; }
echo "[sample-derep] $(grep -c '^>' "$REP_OUT") reps -> $REP_OUT"

# 3) Update global concat (header-dedup) - optimized to avoid intermediate temp file
if [[ -s "$GCONCAT" ]]; then
  # Combine dedup and merge in one step
  awk '
    # First pass: read existing global concat and mark seen headers
    FILENAME==ARGV[1] {
      if (/^>/) {seen[$0]=1; keep=1; print; next}
      if (keep) print
      next
    }
    # Second pass: add new sequences that haven't been seen
    /^>/ {h=$0; if(!(h in seen)){seen[h]=1; keep=1; print; next} keep=0}
    { if(keep) print }
  ' "$GCONCAT" "$REP_OUT" > "$GOUT/.new.concat"
  mv -f "$GOUT/.new.concat" "$GCONCAT"
else
  # First time - just dedup the representative sequences
  awk '
    /^>/ {h=$0; if(!(h in seen)){seen[h]=1; keep=1; print; next} keep=0}
    { if(keep) print }
  ' "$REP_OUT" > "$GCONCAT"
fi

# Optimize sample list update - append and dedup in single operation
{ echo "$SAMPLE"; [[ -f "$GOUT/samples.list" ]] && cat "$GOUT/samples.list"; } | \
  awk '!seen[$0]++' > "$GOUT/.new.samples" && mv -f "$GOUT/.new.samples" "$GOUT/samples.list"
echo "[global] appended $SAMPLE reps to $GCONCAT"

# 4) Trigger global clustering so one run suffices
bash "$ROOT/modules/plasmids/cluster/exec.sh" "$ROOT" _ "${CPUS:-16}" "${FORCE:-0}" || true
