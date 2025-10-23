#!/usr/bin/env bash
# run_plasme_hp.sh
# PLASMe high-precision + Filter: Länge (>= MIN_LEN) UND Ambiguität (Union der amb_region-Intervalle)
# Usage: ./run_plasme_hp.sh ROOT SAMPLE [CPUS] [FORCE]
set -euo pipefail

ROOT="${1:?}"; SAMPLE="${2:?}"; CPUS="${3:-16}"; FORCE="${4:-0}"

# -------- Tunables --------
: "${MIN_LEN:=0}"         # Mindestlänge in bp (empfohlen)
: "${MAX_AMBIG_FRAC:=0.20}"  # max. Anteil ambiger Basen (Union) relativ zur Contig-Länge

# -------- Tools --------
AWK="${AWK:-gawk}"
command -v "$AWK" >/dev/null 2>&1 || { echo "[plasme] gawk not found; set AWK=awk wenn nötig"; exit 2; }
export LC_ALL=C

# -------- Inputs --------
IN_SPADES="$ROOT/SRA/assemblies/spades/$SAMPLE/contigs.fasta"
IN_MPS="$ROOT/SRA/assemblies/metaplasmidspades/$SAMPLE/contigs.fasta"
if   [[ -s "$IN_SPADES" ]]; then IN="$IN_SPADES"
elif [[ -s "$IN_MPS"   ]]; then IN="$IN_MPS"
else  echo "[plasme] no input"; exit 1; fi

DB="$ROOT/refdata/plasme/DB"; [[ -d "$DB" ]] || { echo "[plasme] DB missing: $DB"; exit 2; }

OUT="$ROOT/SRA/plasmids/$SAMPLE/plasme"; mkdir -p "$OUT"
OUT_FASTA_HP="$OUT/plasme.hp.fasta"
FLAG="$OUT/.hp_done"; [[ "$FORCE" == "1" ]] && rm -f "$FLAG"

# -------- Step 1: PLASMe high-precision --------
if [[ ! -f "$FLAG" ]]; then
  echo "[plasme] running PLASMe high-precision..."
  PY="$ROOT/external/PLASMe/PLASMe.py"; command -v python >/dev/null || { echo "[plasme] python missing"; exit 3; }
  [[ -f "$PY" ]] || { echo "[plasme] PLASMe.py not found: $PY"; exit 4; }
  python "$PY" "$IN" "$OUT_FASTA_HP" -d "$DB" -m high-precision -t "$CPUS"
  touch "$FLAG"
else
  echo "[plasme] skip PLASMe (cached). FORCE=1 zum Neu-lauf."
fi

# -------- Outputs --------
REPORT="$(ls -1 "$OUT"/*.tsv "$OUT"/*.csv 2>/dev/null | head -n1 || true)"
PRED_FASTA="$(ls -1 "$OUT"/*.fa "$OUT"/*.fna "$OUT"/*.fasta 2>/dev/null | head -n1 || true)"
[[ -n "${PRED_FASTA:-}" && -s "$PRED_FASTA" ]] || { echo "[plasme] no predicted fasta"; exit 5; }

FILTDIR="$OUT/filtered_hp"; mkdir -p "$FILTDIR"
OUT_FASTA="$FILTDIR/plasme.hp.filtered.fasta"
OUT_TSV="$FILTDIR/plasme.hp.filtered.tsv"
TMP_KEEP="$FILTDIR/keep.ids"; : > "$TMP_KEEP"
DBG="$FILTDIR/debug_parsing.tsv"

# Delimiter auto-detect
if [[ -n "${REPORT:-}" && -s "$REPORT" ]]; then
  header="$(head -n1 "$REPORT")"
  if printf '%s' "$header" | grep -q $'\t'; then SEP=$'\t'
  elif printf '%s' "$header" | grep -q ','; then SEP=',' 
  else SEP=$'\t'; fi
  echo "[plasme] delimiter: $( [[ "$SEP" = $'\t' ]] && echo TAB || echo COMMA )"
fi

# -------- Filter: Länge + Ambiguität (Intervall-Union) --------
if [[ -n "${REPORT:-}" && -s "$REPORT" ]]; then
  "$AWK" -v MINL="$MIN_LEN" -v MAXA="$MAX_AMBIG_FRAC" -v SEP="$SEP" -v DBG="$DBG" '
    # Vereinige a-b,c-d,...; liefere UNION-Länge (inklusive)
    function amb_union_bp(ranges,   n,i,tok,parts,a,b,idx,ns,ne,curS,curE,total,t){
      gsub(/\r/,"",ranges); gsub(/[ \t]/,"",ranges)
      if (ranges=="" || ranges=="-") return 0
      n=split(ranges, tok, /,/); idx=0
      for(i=1;i<=n;i++){
        if (tok[i] ~ /^[0-9]+-[0-9]+$/) {
          split(tok[i], parts, /-/); a=parts[1]+0; b=parts[2]+0
          if (b<a){ t=a; a=b; b=t }
          if (b>=a){ idx++; ns[idx]=a; ne[idx]=b }
        }
      }
      if (idx==0) return 0
      # sort by start
      for(i=2;i<=idx;i++){ a=ns[i]; b=ne[i]; j=i-1
        while(j>=1 && ns[j]>a){ ns[j+1]=ns[j]; ne[j+1]=ne[j]; j-- }
        ns[j+1]=a; ne[j+1]=b
      }
      # merge
      curS=ns[1]; curE=ne[1]; total=0
      for(i=2;i<=idx;i++){
        if (ns[i] <= curE+1) { if (ne[i] > curE) curE = ne[i] }
        else { total += (curE - curS + 1); curS = ns[i]; curE = ne[i] }
      }
      total += (curE - curS + 1)
      return total
    }
    BEGIN{
      FS=SEP; OFS="\t"; IGNORECASE=1
      idc=lenc=ambFracC=ambRegC=0
      print "id","length","amb_bp_union","amb_frac","pass_len","pass_amb","kept" > DBG
    }
    NR==1{
      hdr=0
      for(i=1;i<=NF;i++){
        h=$i; gsub(/\r/,"",h)
        if(h ~ /[A-Za-z]/) hdr=1
        if(h ~ /^(id|contig.*id|sequence.*id)$/) idc=i
        if(h ~ /^(len|length|contig.?len|size|bp)$/) lenc=i
        if((h ~ /(amb|ambig|ambiguous)/) && (h ~ /(frac|ratio|percent|share)/)) ambFracC=i
        if(h ~ /^amb[_ ]?region(s)?$/ || (h ~ /(amb|ambig)/ && h ~ /region/))   ambRegC=i
      }
      if(hdr) next
    }
    {
      id  = (idc?  $idc  : $1)

      # Länge: Spalte, sonst col2, sonst aus ID "_length_####"
      len = -1
      if(lenc)                 len = $lenc + 0
      else if($2 ~ /^[0-9]+$/) len = $2 + 0
      else if(match(id, /_length_([0-9]+)/, m)) len = m[1] + 0

      pass_len = (len < 0 ? 1 : (len >= MINL))

      amb_bp=0; amb=0; pass_amb=1
      if (ambFracC) {
        amb = $ambFracC + 0
        if (amb > 1.000001) amb /= 100.0
        pass_amb = (amb <= MAXA)
      } else if (ambRegC && len > 0) {
        amb_bp = amb_union_bp($ambRegC)
        amb = amb_bp / len
        pass_amb = (amb <= MAXA)
      }

      kept = (pass_len && pass_amb) ? 1 : 0
      print id, len, amb_bp, amb, pass_len, pass_amb, kept >> DBG
      if (kept) print id
    }
  ' "$REPORT" > "$TMP_KEEP"

  [[ -s "$TMP_KEEP" ]] || {
    echo "[plasme] warning: parsing brachte keine IDs → alles behalten"
    grep "^>" "$PRED_FASTA" | sed 's/^>//; s/[ \t].*//' > "$TMP_KEEP"
  }
else
  echo "[plasme] kein Report → alles behalten"
  grep "^>" "$PRED_FASTA" | sed 's/^>//; s/[ \t].*//' > "$TMP_KEEP"
fi

# -------- FASTA filtern --------
"$AWK" '
  BEGIN{ while((getline k<ARGV[1])>0){keep[k]=1}; close(ARGV[1]); ARGV[1]="" }
  /^>/{ id=substr($0,2); sub(/[ \t].*/,"",id); p=(id in keep) }
  { if(p) print }
' "$TMP_KEEP" "$PRED_FASTA" > "$OUT_FASTA"

# -------- Tabelle filtern --------
if [[ -n "${REPORT:-}" && -s "$REPORT" ]]; then
  "$AWK" -v SEP="$SEP" '
    BEGIN{ FS=SEP; OFS="\t"; IGNORECASE=1; idc=0
           while((getline k<ARGV[1])>0){keep[k]=1}; close(ARGV[1]); ARGV[1]="" }
    FNR==1{ for(i=1;i<=NF;i++) if($i ~ /^id$|^contig.*id|^sequence.*id$/) idc=i; print; next }
    { id=(idc?$idc:$1); if(id in keep) print }
  ' "$TMP_KEEP" "$REPORT" > "$OUT_TSV"
fi

# -------- Summary --------
n_all=$(grep -c "^>" "$PRED_FASTA" || true)
n_keep=$(grep -c "^>" "$OUT_FASTA" || true)
echo "[plasme] kept $n_keep / $n_all  ->  $OUT_FASTA"
[[ -s "${OUT_TSV:-/dev/null}" ]] && echo "[plasme] filtered table -> $OUT_TSV"
echo "[plasme] debug table -> $DBG"
