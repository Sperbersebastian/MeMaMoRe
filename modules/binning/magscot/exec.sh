#!/usr/bin/env bash
set -euo pipefail
# env: ROOT SAMPLE
INROOT="$ROOT/SRA/binning"
OUT="$INROOT/magscot/$SAMPLE"
LOG="$ROOT/logs/binning_magscot_${SAMPLE}.log"
mkdir -p "$OUT"
echo "[magscot] start" >"$LOG"

# 0) inputs: collect bins from available binners (MetaBAT2 for now)
BIN_SOURCES=()
for src in metabat2 comebin binny; do
  d="$INROOT/$src/$SAMPLE"
  if compgen -G "$d/bin*.fa*" > /dev/null; then BIN_SOURCES+=("$d"); fi
done
[[ ${#BIN_SOURCES[@]} -gt 0 ]] || { echo "[magscot] no input bins" | tee -a "$LOG"; exit 0; }

# 1) build contigs_to_bin (NO header): BinID \t ContigID \t Tool
CTB="$OUT/contigs_to_bin.tsv"; : > "$CTB"

for tool in metabat2 comebin binny; do
  d="$INROOT/$tool/$SAMPLE"
  if compgen -G "$d/bin*.fa*" >/dev/null; then
    for f in "$d"/bin*.fa*; do
      [[ -s "$f" ]] || continue
      binlabel="$(basename "${f%.*}")"        # e.g., bin.1
      awk -v bin="$binlabel" -v tool="$tool" '
        /^>/ { gsub(/^>/,""); split($0,a," "); print bin "\t" a[1] "\t" tool }
      ' "$f" >> "$CTB"
    done
  fi
done
[[ -s "$CTB" ]] || { echo "[magscot] empty contigs_to_bin" | tee -a "$LOG"; exit 0; }

# 2) build per-sample HMM hits if MAGSCOT_HMM not set (Prodigal + HMMER)
if [[ -z "${MAGSCOT_HMM:-}" ]]; then
  echo "[magscot] building HMM hits (prodigal + hmmsearch)" | tee -a "$LOG"
  ASM="$ROOT/SRA/assemblies/contig_qc/$SAMPLE/contigs.filtered.fasta"
  [[ -s "$ASM" ]] || { echo "[magscot] missing assembly $ASM" | tee -a "$LOG"; exit 0; }
  FAA="$OUT/${SAMPLE}.prodigal.faa"
  micromamba run -n env_mags prodigal -p meta -i "$ASM" -a "$FAA" -o /dev/null >>"$LOG" 2>&1
  [[ -s "$FAA" ]] || { echo "[magscot] prodigal failed" | tee -a "$LOG"; exit 0; }
  # use GTDB r207 marker HMMs shipped with MAGScoT repo
  TIGR="$MAGSCOT_DIR/hmm/gtdbtk_rel207_tigrfam.hmm"
  PFAM="$MAGSCOT_DIR/hmm/gtdbtk_rel207_Pfam-A.hmm"
  [[ -s "$TIGR" && -s "$PFAM" ]] || { echo "[magscot] GTDB r207 HMMs not found in $MAGSCOT_DIR/hmm" | tee -a "$LOG"; exit 0; }
  TOUT="$OUT/${SAMPLE}.hmm.tigr.hit.out"
  POUT="$OUT/${SAMPLE}.hmm.pfam.hit.out"
  micromamba run -n env_mags hmmsearch -o "$OUT/tigr.log"  --tblout "$TOUT"  --noali --notextw --cut_nc --cpu 8 "$TIGR" "$FAA" >>"$LOG" 2>&1
  micromamba run -n env_mags hmmsearch -o "$OUT/pfam.log"  --tblout "$POUT"  --noali --notextw --cut_nc --cpu 8 "$PFAM" "$FAA" >>"$LOG" 2>&1
  # merge to MAGSCOT_HMM format
  T_TAB="$OUT/${SAMPLE}.tigr"; P_TAB="$OUT/${SAMPLE}.pfam"; HMM_MERGED="$OUT/${SAMPLE}.hmm"
  awk '($1 !~ /^#/){print $1"\t"$3"\t"$5}' "$TOUT" > "$T_TAB" || true
  awk '($1 !~ /^#/){print $1"\t"$4"\t"$5}' "$POUT" > "$P_TAB" || true
  cat "$P_TAB" "$T_TAB" > "$HMM_MERGED"
  MAGSCOT_HMM="$HMM_MERGED"
  echo "[magscot] HMM merged -> $MAGSCOT_HMM" | tee -a "$LOG"
fi

# 3) run MAGScoT (inside OUT so outputs land here)
set +e
( cd "$OUT" && micromamba run -n env_mags Rscript "$MAGSCOT_DIR/MAGScoT.R" \
    -i "$CTB" --hmm "$MAGSCOT_HMM" ) >>"$LOG" 2>&1
rc=$?; set -e
if [[ $rc -ne 0 ]]; then echo "[magscot] failed (see $LOG)"; exit 0; fi

# 4) collect refined bins
# collect refined bins
mkdir -p "$OUT/bins"; shopt -s nullglob
mv "$OUT"/*_bin*.fa* "$OUT/bins/" 2>/dev/null || true

if compgen -G "$OUT/bins/*" >/dev/null; then
  echo "[magscot] refined bins -> $OUT/bins"
else
  echo "[magscot] no refined outputs found; falling back to MetaBAT2 bins"
  SRC="$INROOT/metabat2/$SAMPLE"
  if compgen -G "$SRC/bin*.fa*" >/dev/null; then
    cp "$SRC"/bin*.fa* "$OUT/bins/" 2>/dev/null || true
    printf "REFINEMENT_SKIPPED\n" > "$OUT/bins/README.txt"
    echo "[magscot] fallback copied -> $OUT/bins"
  else
    echo "[magscot] no MetaBAT2 bins to fallback to"
  fi
fi
