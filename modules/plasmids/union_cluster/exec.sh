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
  tmp="$UOUT/.gather.$$.fa"; : > "$tmp"

  # viralVerify
  [[ -s "$BASE/viralverify/Prediction_results_fasta/contigs_plasmid.fasta" ]] \
    && awk -v p="vv" '/^>/{sub(/^>/,">"p"|")}1' \
       "$BASE/viralverify/Prediction_results_fasta/contigs_plasmid.fasta" >> "$tmp"
  [[ -s "$BASE/viralverify/Prediction_results_fasta/contigs_plasmid_uncertain.fasta" ]] \
    && awk -v p="vv_uncertain" '/^>/{sub(/^>/,">"p"|")}1' \
       "$BASE/viralverify/Prediction_results_fasta/contigs_plasmid_uncertain.fasta" >> "$tmp"

  # geNomad
  [[ -s "$BASE/genomad/contigs_summary/contigs_plasmid.fna" ]] \
    && awk -v p="genomad" '/^>/{sub(/^>/,">"p"|")}1' \
       "$BASE/genomad/contigs_summary/contigs_plasmid.fna" >> "$tmp"

  # PLASMe (prefer filtered)
  if [[ -s "$BASE/plasme/filtered_hp/plasme.hp.filtered.fasta" ]]; then
    awk -v p="plasme" '/^>/{sub(/^>/,">"p"|")}1' \
      "$BASE/plasme/filtered_hp/plasme.hp.filtered.fasta" >> "$tmp"
  else
    pf="$(ls -1 "$BASE/plasme/"*.fa "$BASE/plasme/"*.fna "$BASE/plasme/"*.fasta 2>/dev/null | head -n1 || true)"
    [[ -n "${pf:-}" ]] && awk -v p="plasme" '/^>/{sub(/^>/,">"p"|")}1' "$pf" >> "$tmp"
  fi

  # MOB-recon
  for f in "$BASE"/mobrecon/*reconstructed*.fasta "$BASE"/mobrecon/*plasmid*.fasta; do
    [[ -s "$f" ]] && awk -v p="mobrecon" '/^>/{sub(/^>/,">"p"|")}1' "$f" >> "$tmp"
  done

  # length filter
  awk -v MIN="$MIN_LEN" '
    /^>/{
      if (seqlen>=MIN && hdr!=""){print hdr; print seq}
      hdr=$0; seq=""; seqlen=0; next
    }
    {seqlen+=length($0); seq=seq $0}
    END{ if (seqlen>=MIN && hdr!=""){print hdr; print seq} }
  ' "$tmp" > "$RAW"
  rm -f "$tmp"

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
  TMP_IDS="$SOUT/rep.ids"
  awk '
    BEGIN{picked=0}
    /^[[:space:]]*$/ {next}
    /^#/ {next}
    /^[>]*[Cc]luster[[:space:]]*[0-9]+/ {picked=0; next}
    { if(!picked){ id=$1; sub(/^[>]/,"",id); print id; picked=1 } }
  ' "$SOUT/clusters.txt" > "$TMP_IDS"

  awk 'BEGIN{
         while((getline k<ARGV[1])>0){want[k]=1}
         close(ARGV[1]); ARGV[1]=""
       }
       /^>/{
         id=substr($0,2); sub(/[ \t].*$/,"",id)
         keep=(id in want)
       }
       { if(keep) print }
  ' "$TMP_IDS" "$LRAW" > "$REP_OUT"
fi
[[ -s "$REP_OUT" ]] || { echo "[err] no representatives"; exit 5; }
echo "[sample-derep] $(grep -c '^>' "$REP_OUT") reps -> $REP_OUT"

# 3) Update global concat (header-dedup)
tmp="$GOUT/.tmp.append.$$.fa"; : > "$tmp"
awk '
  /^>/ {h=$0; if(!(h in seen)){seen[h]=1; keep=1; print; next} keep=0}
  { if(keep) print }
' "$REP_OUT" > "$tmp"

if [[ -s "$GCONCAT" ]]; then
  cat "$GCONCAT" "$tmp" | awk '
    /^>/ {h=$0; if(!(h in seen)){seen[h]=1; keep=1; print; next} keep=0}
    { if(keep) print }
  ' > "$GOUT/.new.concat"
  mv -f "$GOUT/.new.concat" "$GCONCAT"
else
  mv -f "$tmp" "$GCONCAT"
fi
rm -f "$tmp"

echo "$SAMPLE" >> "$GOUT/samples.list"
awk '!seen[$0]++' "$GOUT/samples.list" > "$GOUT/.new.samples" && mv -f "$GOUT/.new.samples" "$GOUT/samples.list"
echo "[global] appended $SAMPLE reps to $GCONCAT"

# 4) Trigger global clustering so one run suffices
bash "$ROOT/modules/plasmids/cluster/exec.sh" "$ROOT" _ "${CPUS:-16}" "${FORCE:-0}" || true
