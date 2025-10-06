#!/usr/bin/env bash
# MAGScoT runner with auto map build, robust parsing, and CTB cleanup
# env: ROOT SAMPLE [CPUS]
set -euo pipefail
: "${ROOT:?}"; : "${SAMPLE:?}"; : "${CPUS:=8}"

OUT="$ROOT/SRA/binning/magscot/$SAMPLE"
mkdir -p "$OUT"

CONTIGS="${CONTIGS:-$ROOT/SRA/assemblies/spades/$SAMPLE/contigs.fasta}"
[[ -s "$CONTIGS" ]] || { echo "[magscot] $SAMPLE no contigs"; exit 0; }

mk_map_from_bins(){ # $1=dir $2=glob $3=tag $4=outfile
  local d="$1" g="$2" tag="$3" out="$4"
  [[ -s "$out" ]] && return 0
  [[ -d "$d" ]] || return 0
  : > "$out"; shopt -s nullglob
  local had=0
  for fa in "$d"/$g; do
    [[ -s "$fa" ]] || continue
    had=1
    bin_base="$(basename "${fa%.*}")"
    bin_id="${tag}|${bin_base}"
    awk -v b="$bin_id" -v t="$tag" '
      BEGIN{RS=">";FS="\n";OFS="\t"}
      NR>1{
        header=$1
        sub(/^[[:space:]]+/,"",header)
        contig=header
        sub(/[[:space:]].*$/,"",contig)
        if(length(contig)>0) print b,contig,t
      }' "$fa" >> "$out"
  done
  (( had )) || : > "$out"
}

# per-tool paths
MB_DIR="$ROOT/SRA/binning/metabat2/$SAMPLE"
BN_DIR="$ROOT/SRA/binning/binny/$SAMPLE/bins"
CB_DIR="$ROOT/SRA/binning/comebin/$SAMPLE/comebin_res/comebin_res_bins"

MB_MAP="$ROOT/SRA/binning/metabat2/$SAMPLE/contigs_to_bin.with_set.tsv"
BN_MAP="$ROOT/SRA/binning/binny/$SAMPLE/contigs_to_bin.with_set.tsv"
CB_MAP="$ROOT/SRA/binning/comebin/$SAMPLE/contigs_to_bin.with_set.tsv"

mk_map_from_bins "$MB_DIR" "*.fa*" "metabat2" "$MB_MAP"
mk_map_from_bins "$BN_DIR" "*.fa*" "binny"    "$BN_MAP"
mk_map_from_bins "$CB_DIR" "*.fa*" "comebin"  "$CB_MAP"

# combine → CTB.raw
RAW="$OUT/contigs_to_bin.raw.tsv"
: > "$RAW"
for f in "$MB_MAP" "$BN_MAP" "$CB_MAP"; do [[ -s "$f" ]] && cat "$f" >> "$RAW"; done
[[ -s "$RAW" ]] || { echo "[magscot] no bin inputs"; exit 0; }

# cleanup:
CTB="$OUT/contigs_to_bin.tsv"
awk -F'\t' 'BEGIN{OFS="\t"} NF==3{k=$3 FS $2; if(!(k in seen)){seen[k]=1; print $0}}' "$RAW" \
| awk -F'\t' 'BEGIN{OFS="\t"} {c[$1]++; rows[NR]=$0} END{for(i=1;i<=NR;i++){split(rows[i],f,"\t"); if(c[f[1]]>=2) print rows[i]}}' \
> "$CTB"

# sanity checks
awk -F'\t' 'NF!=3{bad++} END{if(bad){print "[magscot] bad rows:",bad; exit 2}}' "$CTB"
bins_n=$(cut -f1 "$CTB" | sort -u | wc -l)
contigs_n=$(cut -f2 "$CTB" | sort -u | wc -l)
echo "[magscot] CTB bins=$bins_n contigs=$contigs_n"
if (( contigs_n < bins_n )); then
  echo "[magscot] contigs < bins after cleanup → tighten filter or investigate inputs"; exit 2
fi

# HMM inputs
MAGS_ROOT="${MAGS_ROOT:-$ROOT/tools/MAGScoT}"
HMM_TIGR="$MAGS_ROOT/hmm/gtdbtk_rel207_tigrfam.hmm"
HMM_PFAM="$MAGS_ROOT/hmm/gtdbtk_rel207_Pfam-A.hmm"
SCRIPT="$MAGS_ROOT/MAGScoT.R"
if [[ ! -s "$SCRIPT" || ! -s "$HMM_TIGR" || ! -s "$HMM_PFAM" ]]; then
  echo "[magscot] missing MAGScoT repo or HMMs under $MAGS_ROOT; skipping"; exit 0
fi

TMP="$OUT/tmp"; mkdir -p "$TMP"
FAA="$OUT/prodigal.faa"

# --- Step 4: Prodigal + HMMsearch with skip ---
if [[ ! -s "$FAA" ]]; then
  micromamba run -n env_binning prodigal -i "$CONTIGS" -p meta -a "$FAA" -d "$TMP/prodigal.ffn" -o "$TMP/prodigal.log"
fi

if [[ ! -s "$TMP/tigr.tbl" ]]; then
  micromamba run -n env_binning hmmsearch -o "$TMP/tigr.out" --tblout "$TMP/tigr.tbl" --noali --notextw --cut_nc --cpu "$CPUS" "$HMM_TIGR" "$FAA"
fi
if [[ ! -s "$TMP/pfam.tbl" ]]; then
  micromamba run -n env_binning hmmsearch -o "$TMP/pfam.out" --tblout "$TMP/pfam.tbl" --noali --notextw --cut_nc --cpu "$CPUS" "$HMM_PFAM" "$FAA"
fi

HMM_MAP="$OUT/example.hmm"
if [[ ! -s "$HMM_MAP" ]]; then
  awk 'BEGIN{OFS="\t"} $1!~/^#/{print $4,$1,$14}' "$TMP/tigr.tbl" > "$TMP/tigr.map"
  awk 'BEGIN{OFS="\t"} $1!~/^#/{print $4,$1,$14}' "$TMP/pfam.tbl" > "$TMP/pfam.map"
  cat "$TMP/pfam.map" "$TMP/tigr.map" > "$HMM_MAP"
  cols="$(awk -F'\t' 'NF{print NF; exit}' "$HMM_MAP")"
  [[ "$cols" -eq 3 ]] || { echo "[magscot] bad HMM map ($cols cols) → $HMM_MAP"; exit 1; }
fi

# --- Step 5: Final Rscript with skip ---
if [[ ! -d "$OUT/MAGScoT" || -z "$(ls -A "$OUT/MAGScoT" 2>/dev/null)" ]]; then
  micromamba run -n env_binning Rscript "$SCRIPT" -i "$CTB" --hmm "$HMM_MAP" -o "$OUT/MAGScoT"
else
  echo "[magscot] skip Rscript, results exist in $OUT/MAGScoT"
fi

echo "[magscot] done → $OUT"
