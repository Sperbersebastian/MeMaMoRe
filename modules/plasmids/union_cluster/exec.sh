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

# config (merged PARAMS_YAML from bin/main.sh; defaults when run standalone)
# shellcheck source=bin/lib/cfg.sh
source "$ROOT/bin/lib/cfg.sh"
MIN_LEN="${MIN_LEN:-$(cfg modules.plasmids.min_len 2000)}"
CONSENSUS="$(cfg modules.plasmids.consensus union)"
MIN_OV="$(cfg modules.plasmids.consensus_min_overlap 0.5)"
case "$CONSENSUS" in
  union) MIN_TOOLS=1 ;; min2) MIN_TOOLS=2 ;; min3) MIN_TOOLS=3 ;;
  *) echo "[union] invalid modules.plasmids.consensus '$CONSENSUS' (union|min2|min3)"; exit 2 ;;
esac
STAMP="$UOUT/.consensus"
# rebuild the candidate set when the consensus rule changed since the last run
if [[ -s "$RAW" && "$(cat "$STAMP" 2>/dev/null)" != "$CONSENSUS $MIN_OV" ]]; then
  echo "[union] consensus rule changed -> rebuilding $RAW"; rm -f "$RAW"
fi
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

  # Plasmer (optional 5th predictor; absent when skipped for RAM)
  if [[ -s "$BASE/plasmer/results/${SAMPLE}.plasmer.predPlasmids.fa" ]]; then
    awk -v p="plasmer" '/^>/{sub(/^>/,">"p"|")}1' "$BASE/plasmer/results/${SAMPLE}.plasmer.predPlasmids.fa" >> "$tmp"
  fi

  # MOB-recon
  for f in "$BASE"/mobrecon/*reconstructed*.fasta "$BASE"/mobrecon/*plasmid*.fasta; do
    if [[ -s "$f" ]]; then awk -v p="mobrecon" '/^>/{sub(/^>/,">"p"|")}1' "$f" >> "$tmp"; fi
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

  # consensus: keep candidates confirmed by >= MIN_TOOLS predictors (union = all)
  if [[ "$MIN_TOOLS" -gt 1 && -s "$RAW" ]]; then
    mv -f "$RAW" "$UOUT/plasmids_raw.all.fasta"
    python3 "$ROOT/scripts/plasmid_consensus.py" "$UOUT/plasmids_raw.all.fasta" "$RAW" \
      "$UOUT/plasmid_support.tsv" "$MIN_TOOLS" "$MIN_OV"
  fi
  echo "$CONSENSUS $MIN_OV" > "$STAMP"

  if [[ ! -s "$RAW" ]]; then
    # no plasmid candidates from any tool: normal for some samples
    echo "[union] WARNING: no plasmid sequences >= ${MIN_LEN} bp for $SAMPLE; skipping clustering"
    exit 0
  fi
  echo "[union] $(grep -c '^>' "$RAW") seq -> $RAW"
else
  echo "[union] cached -> $RAW"
fi

# 2) Per-sample clustering (Stampede)
CLSTR="${CLSTR:-$(find "$ROOT/external/Stampede-ClusterGenomes" -type f \( -name 'Cluster_genomes.pl' -o -path '*/Cluster_genomes.pl' \) -print -quit)}"
[[ -n "${CLSTR:-}" && -x "$CLSTR" ]] || { echo "[err] Cluster_genomes.pl not found/executable"; exit 2; }

cp -f "$RAW" "$LRAW"
FLAG="$SOUT/.done.sample.i${ID}.c${COV}"
[[ -f "$FLAG" && "$RAW" -nt "$FLAG" ]] && rm -f "$FLAG"   # candidate set changed
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
