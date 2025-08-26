#!/usr/bin/env bash
set -euo pipefail
# env: ROOT SAMPLE CONTIGS (ignored here), PARAMS_YAML (opt)
OUT="$ROOT/SRA/binning/metabat2/$SAMPLE"
COV="$ROOT/SRA/binning/coverage/$SAMPLE/${SAMPLE}.depth.txt"
ASM_DIR="$ROOT/SRA/assemblies/spades/$SAMPLE"
LOG="$ROOT/logs/binning_metabat2_${SAMPLE}.log"
mkdir -p "$OUT"

# pick the SPAdes assembly used for mapping
RAW_CONTIG="$ASM_DIR/contigs.len1000.fasta"
[[ -s "$RAW_CONTIG" ]] || RAW_CONTIG="$ASM_DIR/contigs.fasta"
[[ -s "$RAW_CONTIG" ]] || { echo "[metabat2] no contigs for $SAMPLE"; exit 0; }

[[ -s "$COV" ]] || { echo "[metabat2] coverage missing for $SAMPLE"; exit 2; }

# threads
THREADS=$(micromamba run -n env_sra_tools python "$ROOT/scripts/merge_params.py" \
  --defaults "$PARAMS_YAML" --get global.threads 2>/dev/null || echo 4)

# Filter contigs >=1500bp (MetaBAT2 hard floor)
FILT="$OUT/contigs.min1500.fasta"
micromamba run -n env_binning seqkit seq -m 1500 "$RAW_CONTIG" > "$FILT"

if [[ ! -s "$FILT" ]]; then
  echo "[metabat2] no contigs >=1500 for $SAMPLE" | tee -a "$LOG"; exit 0;
fi

# filter depth rows to same contigs AND reorder to match FASTA order
KEEP="$OUT/contigs.keep.txt"; grep '^>' "$FILT" | sed 's/^>//; s/ .*//' > "$KEEP"
DEP_FILT="$OUT/${SAMPLE}.depth.min1500.txt"
awk 'NR==1{print; next} FNR==NR{keep[$1]=1; next} ($1 in keep)' "$KEEP" "$COV" > "$DEP_FILT"

# Build JGI depth in EXACT FASTA order; zero-fill missing contigs
DEP_SORT="$OUT/${SAMPLE}.depth.min1500.ordered.txt"
micromamba run -n env_sra_tools python - <<'PY' "$FILT" "$COV" "$DEP_SORT"
import sys, csv
fa, dep_in, dep_out = sys.argv[1:]
# 1) FASTA order
order=[]
with open(fa) as f:
    for line in f:
        if line.startswith('>'):
            order.append(line[1:].split()[0])
# 2) read JGI depth to dict: name -> row (keep header)
rows={}
with open(dep_in, newline='') as f:
    r=csv.reader(f, delimiter='\t')
    header=next(r)
    # detect columns (name, len, mean, var) by position
    for row in r:
        if not row: continue
        rows[row[0]]=row
# 3) write header, then rows in FASTA order; zero-fill if missing
with open(dep_out,'w',newline='') as o:
    w=csv.writer(o, delimiter='\t')
    w.writerow(header)
    # find indices or fall back to positions
    try:
        i_len = header.index('contigLen')
    except ValueError:
        i_len = 1
    try:
        i_mean = header.index('totalAvgDepth')  # jgi often uses totalAvgDepth
    except ValueError:
        i_mean = 2
    try:
        i_var = next(i for i,h in enumerate(header) if h.endswith('-var'))
    except StopIteration:
        i_var = i_mean+1
    for name in order:
        if name in rows:
            w.writerow(rows[name])
        else:
            # zero row: we don't know contigLen from JGI; set 0 (MetaBAT ignores it)
            z = [None]*len(header)
            z[0] = name
            z[i_len]  = '0'
            z[i_mean] = '0.000000'
            z[i_var]  = '0.000000'
            w.writerow(z)
PY

[[ -s "$DEP_SORT" ]] || { echo "[metabat2] ordered depth empty"; exit 0; }


# Run MetaBAT2
micromamba run -n env_binning metabat2 \
  -i "$FILT" \
  -a "$DEP_SORT" \
  -o "$OUT/bin" \
  -t "$THREADS" \
  --minContig 1500 \
  >"$LOG" 2>&1 || true

echo "[metabat2] bins -> $OUT/"
