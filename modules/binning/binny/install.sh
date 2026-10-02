#!/usr/bin/env bash
set -euo pipefail

# ROOT is the project root (defaults to the repo containing this script)
ROOT="${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
BINNY_DIR="$ROOT/external/binny"

TF="$BINNY_DIR/database/hmms/checkm_tf"
PF="$BINNY_DIR/database/hmms/checkm_pf"

# --- prepare support dirs to vendor tiny, persistent assets
SUP_RES="$BINNY_DIR/database/support/Resources/NCBI"
SUP_NLTK="$BINNY_DIR/database/support/nltk_data"
mkdir -p "$SUP_RES" "$SUP_NLTK"

# --- vendor NCBI genetic code file (tiny)
if [[ ! -s "$SUP_RES/gc.prt.dmp" ]]; then
  tmpd="$(mktemp -d)"; pushd "$tmpd" >/dev/null
  curl -fsSL https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz | tar -xz gc.prt
  install -m0644 gc.prt "$SUP_RES/gc.prt.dmp"
  popd >/dev/null; rm -rf "$tmpd"
fi

# --- vendor NLTK models into repo (tiny)
python - <<PY
import nltk, os, sys
base = os.environ.get("SUP_NLTK","$SUP_NLTK")
os.makedirs(base, exist_ok=True)
def have(path): 
    import nltk.data
    try: nltk.data.find(path); return True
    except LookupError: return False
todo=[]
if not have("taggers/averaged_perceptron_tagger_eng/"): todo.append("averaged_perceptron_tagger_eng")
if not have("tokenizers/punkt/"): todo.append("punkt")
for pkg in todo: nltk.download(pkg, download_dir=base, quiet=True)
print("NLTK ready at:", base)
PY

# --- HMM housekeeping: hmmpress + metadata.tsv (idempotent)
prep_hmm_dir() {
  local D="$1"; local HMM=("$D"/checkm_filtered_*.hmm)
  [[ -s "${HMM[0]}" ]] || { echo "[install] missing HMM in $D"; exit 2; }
  # build indices if any are missing
  if ! ls "$D"/checkm_filtered_*.hmm.h3{f,i,m,p} >/dev/null 2>&1; then
    hmmpress "${HMM[0]}"
  fi
  # metadata.tsv in MANTIS format (REF, dummy '|', typed links)
  if [[ ! -s "$D/metadata.tsv" ]]; then
    awk -v OFS='\t' '
      /^NAME[ \t]/{n=$2}
      /^ACC[ \t]/{a=$2}
      /^LENG[ \t]/{len=$2; db=(a ~ /^PF/ ? "pfam" : (a ~ /^TIGR/ ? "tigrfam" : "custom"));
                  if(n!=""){print n, "|", "acc:"a, "db:"db, "length:"len; n=a=len=""}}
    ' "${HMM[0]}" > "$D/metadata.tsv"
  fi
}
prep_hmm_dir "$TF"
prep_hmm_dir "$PF"

echo "[install] binny assets ready:"
echo " - $SUP_RES/gc.prt.dmp"
echo " - $SUP_NLTK (nltk_data)"
echo " - $TF/metadata.tsv + *.h3*"
echo " - $PF/metadata.tsv + *.h3*"
