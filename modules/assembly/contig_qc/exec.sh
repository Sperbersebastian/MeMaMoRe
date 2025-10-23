#!/usr/bin/env bash
set -euo pipefail
# env: ROOT PARAMS_YAML SAMPLE MODE R1 R2 FORCE

# --- ensure dirs ---
LOGS="$ROOT/logs"
EVENTS="$ROOT/api/status.jsonl"
mkdir -p "$LOGS" "$(dirname "$EVENTS")"

# --- params ---
THREADS=$(micromamba run -n env_sra_tools python "$ROOT/scripts/merge_params.py" --defaults "$PARAMS_YAML" --get global.threads)
BREADTH_MIN=$(micromamba run -n env_sra_tools python "$ROOT/scripts/merge_params.py" --defaults "$PARAMS_YAML" --get modules.assembly.qc.breadth_min)
DEPTH_MIN=$(micromamba run -n env_sra_tools python "$ROOT/scripts/merge_params.py" --defaults "$PARAMS_YAML" --get modules.assembly.qc.depth_min)

# --- paths ---
OUTBASE="$ROOT/SRA/assemblies/spades/$SAMPLE"
MAPDIR="$OUTBASE/map"; mkdir -p "$MAPDIR"

QCBASE="$ROOT/SRA/assemblies/contig_qc/$SAMPLE"
mkdir -p "$QCBASE"

keep_tsv="$QCBASE/contigs_kept.tsv"
drop_tsv="$QCBASE/contigs_dropped.tsv"
filtered_fa="$QCBASE/contigs.filtered.fasta"
mapcounts="$QCBASE/mapping_counts.tsv"

ts(){ date -u +"%Y-%m-%dT%H:%M:%SZ"; }
echo "{\"ts\":\"$(ts)\",\"module\":\"assembly\",\"program\":\"contig_qc\",\"sample\":\"$SAMPLE\",\"phase\":\"start\"}" >> "$EVENTS"

# --- assembly fasta (prefer len-filtered) ---
REF="$OUTBASE/contigs.len1000.fasta"; [[ -s "$REF" ]] || REF="$OUTBASE/contigs.fasta"
[[ -s "$REF" ]] || { echo "[contig_qc] no assembly for $SAMPLE" >&2; exit 2; }

# --- choose reads (test fallback to raw if missing) ---
if [[ "${MODE:-full}" == "test" ]]; then
  if [[ -n "${R2:-}" ]]; then
    IN1="$ROOT/SRA/testsets/$SAMPLE/${SAMPLE}_1.10p.fastq.gz"
    IN2="$ROOT/SRA/testsets/$SAMPLE/${SAMPLE}_2.10p.fastq.gz"
  else
    IN1="$ROOT/SRA/testsets/$SAMPLE/${SAMPLE}.10p.fastq.gz"; IN2=""
  fi
  [[ -s "$IN1" ]] || { IN1="$R1"; IN2="${R2:-}"; }
else
  T1="$ROOT/SRA/qc/fastp/$SAMPLE/${SAMPLE}_trimmed_1.fastq.gz"
  T2="$ROOT/SRA/qc/fastp/$SAMPLE/${SAMPLE}_trimmed_2.fastq.gz"
  if [[ -s "$T1" ]]; then
    IN1="$T1"; IN2=""; [[ -s "$T2" ]] && IN2="$T2"
  else
    IN1="$R1"; IN2="${R2:-}"
  fi
fi

# --- map reads (BWA-MEM2) ---
BAM="$MAPDIR/${SAMPLE}.sorted.bam"
if [[ -s "$BAM" && "${FORCE:-0}" != "1" ]]; then
  echo "[contig_qc] Reusing existing BAM for $SAMPLE" >> "$LOGS/map_${SAMPLE}.log"
else
  [[ "${FORCE:-0}" == "1" ]] && rm -f "$MAPDIR"/* || true
  # index reference if needed
  if [[ ! -s "${REF}.0123" && ! -s "${REF}.bwt.2bit.64" ]]; then
    micromamba run -n env_mapping_coverm bwa-mem2 index "$REF" >"$LOGS/map_${SAMPLE}.log" 2>&1
  fi
  if [[ -n "$IN2" ]]; then
    micromamba run -n env_mapping_coverm bash -lc \
      "bwa-mem2 mem -t $THREADS '$REF' '$IN1' '$IN2' | samtools sort -@ $THREADS -o '$BAM'" >>"$LOGS/map_${SAMPLE}.log" 2>&1
  else
    micromamba run -n env_mapping_coverm bash -lc \
      "bwa-mem2 mem -t $THREADS '$REF' '$IN1' | samtools sort -@ $THREADS -o '$BAM'" >>"$LOGS/map_${SAMPLE}.log" 2>&1
  fi
  micromamba run -n env_mapping_coverm samtools index "$BAM" >>"$LOGS/map_${SAMPLE}.log" 2>&1
fi

# --- CoverM per-contig ---
COVALL="$MAPDIR/coverm_contigs.tsv"
if [[ -n "$IN2" ]]; then
  micromamba run -n env_mapping_coverm coverm contig \
    --reference "$REF" -1 "$IN1" -2 "$IN2" \
    --methods covered_fraction mean variance length \
    --min-read-aligned-length 50 --threads "$THREADS" \
    > "$COVALL" 2>>"$LOGS/map_${SAMPLE}.log"
else
  micromamba run -n env_mapping_coverm coverm contig \
    --reference "$REF" --single "$IN1" \
    --methods covered_fraction mean variance length \
    --min-read-aligned-length 50 --threads "$THREADS" \
    > "$COVALL" 2>>"$LOGS/map_${SAMPLE}.log"
fi

# --- filter contigs + write filtered FASTA ---
export REF COVALL keep_tsv drop_tsv filtered_fa BREADTH_MIN DEPTH_MIN
micromamba run -n env_sra_tools python - <<'PY'
import csv, os, math, sys

ref         = os.environ["REF"]
cov         = os.environ["COVALL"]
keep_tsv    = os.environ["keep_tsv"]
drop_tsv    = os.environ["drop_tsv"]
filt_fa     = os.environ["filtered_fa"]
breadth_min = float(os.environ.get("BREADTH_MIN","0.8"))
depth_min   = float(os.environ.get("DEPTH_MIN","2.0"))

if not os.path.exists(cov) or os.path.getsize(cov)==0:
    raise SystemExit(f"[contig_qc] CoverM table empty/missing: {cov}")

with open(cov, newline='') as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)

def find_suffix(*suffixes):
    lower = [h.lower() for h in header]
    for suf in suffixes:
        suf = suf.lower()
        for i,h in enumerate(lower):
            if h.endswith(suf):
                return header[i]
    raise KeyError(f"no column endswith {suffixes}; header={header}")

C_CONTIG  = next((h for h in header if h.lower()=="contig"), None)
if C_CONTIG is None:
    raise KeyError(f"No 'Contig' column in header={header}")

C_BREADTH = find_suffix(" covered fraction")
C_MEAN    = find_suffix(" mean")
C_VAR     = find_suffix(" variance")
C_LEN     = find_suffix(" length")

rows=[]
with open(cov, newline='') as f:
    r=csv.reader(f, delimiter='\t')
    next(r)
    for row in r:
        if not row: continue
        try:
            cols = {name:i for i,name in enumerate(header)}
            name = row[cols[C_CONTIG]].split()[0]
            br   = float(row[cols[C_BREADTH]])
            md   = float(row[cols[C_MEAN]])
            var  = float(row[cols[C_VAR]])
            ln   = int(float(row[cols[C_LEN]]))
            cv   = (math.sqrt(var)/md) if md>0 else float('inf')
        except Exception as e:
            print(f"[warn] skip row: {row} ({e})", file=sys.stderr)
            continue
        reason=[]
        if br < breadth_min: reason.append("breadth")
        if md < depth_min:   reason.append("depth")
        rows.append(dict(contig=name,length=ln,breadth=br,mean=md,variance=var,cv=cv,reason=",".join(reason)))

keep=[r for r in rows if r['reason']==""]
drop=[r for r in rows if r['reason']!=""]

def wtsv(path, arr):
    with open(path, "w", newline='') as o:
        w=csv.writer(o, delimiter='\t')
        w.writerow(["contig","length","breadth","mean","variance","cv","reason"])
        for r in arr:
            w.writerow([r['contig'], r['length'],
                        f"{r['breadth']:.6f}", f"{r['mean']:.6f}",
                        f"{r['variance']:.6f}", f"{r['cv']:.6f}", r['reason']])
wtsv(keep_tsv, keep); wtsv(drop_tsv, drop)

keepl=set(r['contig'] for r in keep)

def stream_fa(path):
    name=None; buf=[]
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None: yield name, ''.join(buf)
                name=line[1:].strip().split()[0]; buf=[]
            else:
                buf.append(line.strip())
        if name is not None: yield name, ''.join(buf)

with open(filt_fa,'w') as o:
    for name,seq in stream_fa(ref):
        if name in keepl:
            o.write(f'>{name}\n')
            for i in range(0,len(seq),80):
                o.write(seq[i:i+80]+'\n')
PY

# --- mapping summaries ---
IDX="$MAPDIR/idxstats.tsv"
micromamba run -n env_mapping_coverm samtools idxstats "$BAM" > "$IDX"
micromamba run -n env_mapping_coverm samtools view -c -F 256 -F 2048 "$BAM" > "$MAPDIR/total_primary.txt"

export mapcounts MAPDIR keep_tsv SAMPLE
micromamba run -n env_sra_tools python - <<'PY'
import csv, os
idx = os.path.join(os.environ["MAPDIR"], "idxstats.tsv")
totp = os.path.join(os.environ["MAPDIR"], "total_primary.txt")
keep_tsv = os.environ["keep_tsv"]
out_tsv  = os.environ["mapcounts"]

kept=set()
with open(keep_tsv) as f:
    r=csv.DictReader(f, delimiter='\t')
    for row in r: kept.add(row['contig'])

total_primary = int(open(totp).read().strip())
mapped_unfiltered = 0; mapped_filtered = 0
with open(idx) as f:
    for line in f:
        ref, length, mapped, unmapped = line.rstrip().split('\t')
        if ref == '*': continue
        m = int(mapped); mapped_unfiltered += m
        if ref in kept: mapped_filtered += m

def prop(a,b): return (a/b) if b>0 else 0.0
with open(out_tsv, "w") as o:
    w=csv.writer(o, delimiter='\t')
    w.writerow(["sample","total_primary","mapped_unfiltered","mapped_filtered","prop_unfiltered","prop_filtered"])
    w.writerow([os.environ["SAMPLE"], total_primary, mapped_unfiltered, mapped_filtered,
                f"{prop(mapped_unfiltered,total_primary):.6f}", f"{prop(mapped_filtered,total_primary):.6f}"])
PY

echo "{\"ts\":\"$(ts)\",\"module\":\"assembly\",\"program\":\"contig_qc\",\"sample\":\"$SAMPLE\",\"phase\":\"done\",\"out\":{\"bam\":\"$BAM\",\"cov\":\"$COVALL\",\"kept\":\"$keep_tsv\",\"drop\":\"$drop_tsv\",\"fasta\":\"$filtered_fa\",\"mapcounts\":\"$mapcounts\"}}" >> "$EVENTS"
