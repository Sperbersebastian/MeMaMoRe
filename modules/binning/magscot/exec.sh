#!/usr/bin/env bash
# MAGScoT runner: nutzt nur die normalisierte RAW-Tabelle (3 Spalten) als Input
# env: ROOT SAMPLE [CPUS]
set -euo pipefail
: "${ROOT:?}"; : "${SAMPLE:?}"; : "${CPUS:=8}"

OUT="$ROOT/SRA/binning/magscot/$SAMPLE"
mkdir -p "$OUT"
MAPDIR="$OUT/maps"; mkdir -p "$MAPDIR"

CONTIGS="${CONTIGS:-$ROOT/SRA/assemblies/spades/$SAMPLE/contigs.fasta}"
[[ -s "$CONTIGS" ]] || { echo "[magscot] $SAMPLE no contigs"; exit 0; }

# ---------- Map-Erzeugung nur falls NICHT vorhanden ----------
mk_map_from_bins(){ # $1=dir $2=glob $3=tag $4=outfile
  local d="$1" g="$2" tag="$3" out="$4"
  [[ -s "$out" ]] && return 0
  [[ -d "$d" ]] || return 0
  : > "$out"
  shopt -s nullglob
  local had=0
  for fa in "$d"/$g; do
    [[ -s "$fa" ]] || continue
    had=1
    local bin_base reader
    bin_base="$(basename "${fa%.*}")"
    case "$fa" in *.gz) reader="zcat";; *) reader="cat";; esac
    $reader "$fa" | awk -v b="${tag}|${bin_base}" -v t="$tag" '
      BEGIN{RS=">";FS="\n";OFS="\t"}
      NR>1{
        header=$1; sub(/^[[:space:]]+/,"",header)
        contig=header; sub(/[[:space:]].*$/,"",contig)
        if(length(contig)>0) print b,contig,t
      }' >> "$out"
  done
  shopt -u nullglob
  (( had )) || : > "$out"
}

# ---------- Modul-Pfade ----------
MB_DIR="$ROOT/SRA/binning/metabat2/$SAMPLE"
BN_DIR="$ROOT/SRA/binning/binny/$SAMPLE/bins"
CB_DIR="$ROOT/SRA/binning/comebin/$SAMPLE/comebin_res/comebin_res_bins"
SB_DIR="$ROOT/SRA/binning/semibin2/$SAMPLE/output_bins"

MB_EXIST_1="$ROOT/SRA/binning/metabat2/$SAMPLE/contigs_to_bin.with_set.tsv"
MB_EXIST_2="$ROOT/SRA/binning/metabat2/$SAMPLE/contigs_to_bin.tsv"
BN_EXIST_1="$ROOT/SRA/binning/binny/$SAMPLE/contigs_to_bin.with_set.tsv"
BN_EXIST_2="$ROOT/SRA/binning/binny/$SAMPLE/contigs_to_bin.tsv"
CB_EXIST_1="$ROOT/SRA/binning/comebin/$SAMPLE/contigs_to_bin.with_set.tsv"
CB_EXIST_2="$ROOT/SRA/binning/comebin/$SAMPLE/contigs_to_bin.tsv"
SB_EXIST_1="$ROOT/SRA/binning/semibin2/$SAMPLE/contigs_to_bin.with_set.tsv"

# Zielkarten nur unter $OUT/maps
MB_MAP="$MAPDIR/metabat2.with_set.tsv"
BN_MAP="$MAPDIR/binny.with_set.tsv"
CB_MAP="$MAPDIR/comebin.with_set.tsv"
SB_MAP="$MAPDIR/semibin2.with_set.tsv"

# Symlink oder on-the-fly Build (Tool-Ordner bleiben unangetastet)
if   [[ -s "$MB_EXIST_1" ]]; then ln -sf "$MB_EXIST_1" "$MB_MAP"
elif [[ -s "$MB_EXIST_2" ]]; then ln -sf "$MB_EXIST_2" "$MB_MAP"
else mk_map_from_bins "$MB_DIR" "bin.*.fa" "metabat2" "$MB_MAP"; fi

if   [[ -s "$BN_EXIST_1" ]]; then ln -sf "$BN_EXIST_1" "$BN_MAP"
elif [[ -s "$BN_EXIST_2" ]]; then ln -sf "$BN_EXIST_2" "$BN_MAP"
else mk_map_from_bins "$BN_DIR" "*.fa*" "binny" "$BN_MAP"; fi

if   [[ -s "$CB_EXIST_1" ]]; then ln -sf "$CB_EXIST_1" "$CB_MAP"
elif [[ -s "$CB_EXIST_2" ]]; then ln -sf "$CB_EXIST_2" "$CB_MAP"
else mk_map_from_bins "$CB_DIR" "*.fa*" "comebin" "$CB_MAP"; fi

if   [[ -s "$SB_EXIST_1" ]]; then ln -sf "$SB_EXIST_1" "$SB_MAP"
else mk_map_from_bins "$SB_DIR" "*.fa" "semibin2" "$SB_MAP"; fi

# ---------- Normalisierung auf 3 Spalten + contig/bin-Korrektur ----------
normalize_map(){ # $1=in $2=tag
  local in="$1" tag="$2" tmp="${1}.norm"
  [[ -s "$in" ]] || { : > "$tmp"; mv -f "$tmp" "$in"; return 0; }
  awk -v t="$tag" -F'\t' '
    BEGIN{OFS="\t"}
    function iscontig(x){
      return (x ~ /^(NODE_|contig|scaffold|tig|chr|NC_|NZ_|gi\||[A-Z]{2}_[A-Z0-9]+|J[A-Z]|C[A-Z]|M[A-Z]|U[A-Z]|K[A-Z]|L[A-Z]|A[A-Z])/)
    }
    /^#/ {next}
    NF>=2{
      a=$1; b=$2; set=(NF>=3 && $3!="")?$3:t
      # ensure a=bin, b=contig
      if (iscontig(a) && !iscontig(b)) { tmp=a; a=b; b=tmp }
      if (iscontig(a) && iscontig(b)) next
      if (index(a, t "|") != 1) a = t "|" a
      print a, b, set
    }' "$in" > "$tmp" && mv -f "$tmp" "$in"
}

normalize_map "$MB_MAP" metabat2
normalize_map "$BN_MAP" binny
normalize_map "$CB_MAP" comebin
normalize_map "$SB_MAP" semibin2

# ---------- RAW normalisiert aufbauen (alle Tools), nur identische Zeilen entfernen ----------
RAW_NORM="$OUT/contigs_to_bin.raw.norm.tsv"
: > "$RAW_NORM"
for f in "$MB_MAP" "$BN_MAP" "$CB_MAP" "$SB_MAP"; do
  if [[ -s "$f" ]]; then cat "$f"; fi
done \
| awk -F'\t' 'BEGIN{OFS="\t"} NF==3{k=$1 FS $2 FS $3; if(!(k in seen)){seen[k]=1; print}}' \
> "$RAW_NORM"

[[ -s "$RAW_NORM" ]] || { echo "[magscot] no bin inputs"; exit 0; }

# ---------- Aufräumen: keine anderen Tabellen behalten ----------
rm -f "$OUT/contigs_to_bin.raw.tsv" "$OUT/contigs_to_bin.tsv" 2>/dev/null || true

# ---------- HMM + MAGScoT ----------
MAGS_ROOT="${MAGS_ROOT:-$ROOT/external/MAGScoT}"
HMM_TIGR="$MAGS_ROOT/hmm/gtdbtk_rel207_tigrfam.hmm"
HMM_PFAM="$MAGS_ROOT/hmm/gtdbtk_rel207_Pfam-A.hmm"
SCRIPT="$MAGS_ROOT/MAGScoT.R"
if [[ ! -s "$SCRIPT" || ! -s "$HMM_TIGR" || ! -s "$HMM_PFAM" ]]; then
  echo "[magscot] missing MAGScoT/HMMs under $MAGS_ROOT; skipping"; exit 0
fi

TMP="$OUT/tmp"; mkdir -p "$TMP"
# central gene calls (modules/assembly/genes); own Prodigal run only as fallback
FAA="$ROOT/SRA/assemblies/genes/$SAMPLE/genes.faa"
[[ -s "$FAA" ]] || FAA="$OUT/prodigal.faa"

[[ -s "$FAA" ]] || micromamba run -n env_binning prodigal -i "$CONTIGS" -p meta -a "$FAA" -d "$TMP/prodigal.ffn" -o "$TMP/prodigal.log"
[[ -s "$TMP/tigr.tbl" && "$TMP/tigr.tbl" -nt "$FAA" ]] || micromamba run -n env_binning hmmsearch -o "$TMP/tigr.out" --tblout "$TMP/tigr.tbl" --noali --notextw --cut_nc --cpu "$CPUS" "$HMM_TIGR" "$FAA"
[[ -s "$TMP/pfam.tbl" && "$TMP/pfam.tbl" -nt "$FAA" ]] || micromamba run -n env_binning hmmsearch -o "$TMP/pfam.out" --tblout "$TMP/pfam.tbl" --noali --notextw --cut_nc --cpu "$CPUS" "$HMM_PFAM" "$FAA"

# MAGScoT --hmm: gene id, marker accession, e-value (tblout: $1 target, $4 query acc, $5 e-value).
# Older runs wrote "marker, gene, 0", which made MAGScoT score no bin at all: rebuild those.
HMM_MAP="$OUT/example.hmm"
if [[ -s "$HMM_MAP" ]] && awk -F'\t' 'NR==1{exit !($1 ~ /^(PF|TIGR)[0-9]/)}' "$HMM_MAP"; then
  echo "[magscot] rebuilding $HMM_MAP (old column order) and discarding MAGScoT results built from it"
  rm -f "$HMM_MAP" "$OUT"/MAGScoT.*
fi
[[ -s "$HMM_MAP" && "$HMM_MAP" -nt "$TMP/tigr.tbl" ]] || rm -f "$HMM_MAP" "$OUT"/MAGScoT.*
if [[ ! -s "$HMM_MAP" ]]; then
  awk 'BEGIN{OFS="\t"} $1!~/^#/{print $1,$4,$5}' "$TMP/tigr.tbl" > "$TMP/tigr.map"
  awk 'BEGIN{OFS="\t"} $1!~/^#/{print $1,$4,$5}' "$TMP/pfam.tbl" > "$TMP/pfam.map"
  cat "$TMP/pfam.map" "$TMP/tigr.map" > "$HMM_MAP"
  cols="$(awk -F'\t' 'NF{print NF; exit}' "$HMM_MAP")"
  [[ "$cols" -eq 3 ]] || { echo "[magscot] bad HMM map ($cols cols) → $HMM_MAP"; exit 1; }
fi

INPUT_FOR_MAGS="$RAW_NORM"
# MAGScoT writes <prefix>.refined.contig_to_bin.out (binnew \t contig), <prefix>.refined.out, <prefix>.scores.out
REFINED="$OUT/MAGScoT.refined.contig_to_bin.out"
# results are only reused for the same candidate bins (e.g. not after a new binner)
IN_SUM="$(md5sum < "$INPUT_FOR_MAGS" | cut -d' ' -f1)"
if [[ -s "$REFINED" && "$(cat "$OUT/.magscot_input.md5" 2>/dev/null)" != "$IN_SUM" ]]; then
  echo "[magscot] candidate bins changed -> recomputing"
  rm -f "$OUT"/MAGScoT.*
fi
if [[ ! -s "$REFINED" ]]; then
  # MAGScoT's bin-merging step crashes when no bin pair qualifies for merging
  # (left_join on an empty table); retry without merging in that case.
  if ! micromamba run -n env_binning Rscript "$SCRIPT" -i "$INPUT_FOR_MAGS" --hmm "$HMM_MAP" -o "$OUT/MAGScoT"; then
    echo "[magscot] MAGScoT failed, retrying with --skip_merge_bins"
    micromamba run -n env_binning Rscript "$SCRIPT" -i "$INPUT_FOR_MAGS" --hmm "$HMM_MAP" -o "$OUT/MAGScoT" --skip_merge_bins
  fi
else
  echo "[magscot] skip Rscript, results exist: $REFINED"
fi
# With bins from a single binning tool MAGScoT only scores them ("Scoring bins
# only") and writes no refined table: keep that tool's bins whose MAGScoT
# score passes the same threshold as refinement (MAGSCOT_MIN_SCORE, 0.5).
SCORES="$OUT/MAGScoT.scores.out"
if [[ ! -s "$REFINED" && -s "$SCORES" ]]; then
  awk -F'\t' -v min="${MAGSCOT_MIN_SCORE:-0.5}" '
    FNR==NR { if (FNR==1) { for (i=1;i<=NF;i++) if ($i=="score") sc=i; next }
              if (sc && $sc+0 >= min) keep[$1]=1; next }
    FNR==1  { print "binnew\tcontig" }
    ($1 in keep) { b=$1; gsub(/[|\/]/, "_", b); print b "\t" $2 }
  ' "$SCORES" "$RAW_NORM" > "$REFINED.tmp"
  if [[ $(wc -l < "$REFINED.tmp") -gt 1 ]]; then
    mv -f "$REFINED.tmp" "$REFINED"
    echo "[magscot] single binning tool: kept bins with score >= ${MAGSCOT_MIN_SCORE:-0.5} as refined set"
  else
    rm -f "$REFINED.tmp"
  fi
fi
echo "$IN_SUM" > "$OUT/.magscot_input.md5"
[[ -s "$REFINED" ]] || { echo "[magscot] no refined bins ($REFINED missing)"; exit 0; }

# ---------- refined bins -> bins/<bin>.fa (input for CheckM2, barrnap/tRNAscan, GTDB-Tk) ----------
# MAGScoT prefixes refined bin names with the -o value, i.e. the full output path: keep the basename
BINS_DIR="$OUT/bins"
rm -rf "$BINS_DIR"; mkdir -p "$BINS_DIR"
awk -v dir="$BINS_DIR" '
  FNR==NR { if (FNR>1 && NF>=2) { b=$1; sub(/.*\//, "", b); bin[$2]=b }; next }
  /^>/    { id=substr($1,2); out=(id in bin) ? dir "/" bin[id] ".fa" : ""; if (out!="") print > out; next }
  out!="" { print > out }
' "$REFINED" "$CONTIGS"
# refined contig->bin table for the GUI
awk 'BEGIN{OFS="\t"; print "contig","bin"} NR>1 && NF>=2 {b=$1; sub(/.*\//, "", b); print $2,b}' "$REFINED" > "$OUT/contigs_to_bin.tsv"
echo "[magscot] $(find "$BINS_DIR" -name '*.fa' | wc -l) refined bins -> $BINS_DIR"

echo "[magscot] done → $OUT"
