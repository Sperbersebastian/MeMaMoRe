#!/usr/bin/env bash
# End-to-end test on a simulated community with known composition.
#
#   tests/e2e/run_e2e.sh [--config config/sim.yaml] [--cpus N] [--skip-binning] [--keep]
#
# Steps: simulate reads (modules/sim) -> ingest -> qc -> assembly -> [binning] -> plasmids
#        -> align contigs and plasmid predictions to the known references -> metrics.
# Results: tests/e2e/results/<SET>/{metrics.tsv,metrics.json,versions.tsv,timings.tsv}
#
# Needs micromamba plus the databases the modules use (geNomad, PLASMe, CheckM2, GTDB-Tk
# for binning). Not run in CI: it downloads references and takes hours on a workstation.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CFG="$ROOT/config/sim.yaml"
CPUS="${CPUS:-8}"
SKIP_BINNING=0
KEEP=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CFG="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"; shift 2 ;;
    --cpus) CPUS="$2"; shift 2 ;;
    --skip-binning) SKIP_BINNING=1; shift ;;
    --keep) KEEP=1; shift ;;
    -h|--help) sed -n 2,13p "$0"; exit 0 ;;
    *) echo "[e2e] unknown option: $1" >&2; exit 2 ;;
  esac
done
export CFG CPUS

MAIN="$ROOT/bin/main.sh"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
command -v micromamba >/dev/null || { echo "[e2e] micromamba not found" >&2; exit 127; }

# set name from the config (plain grep: yq lives inside sim_env, which may not exist yet)
SET="$(grep -E '^set_name:' "$CFG" | head -n1 | sed -E 's/^set_name:[[:space:]]*"?([^"#]*)"?.*/\1/' | xargs)"
[[ -n "$SET" ]] || { echo "[e2e] set_name missing in $CFG" >&2; exit 2; }

OUT="$ROOT/tests/e2e/results/$SET"; mkdir -p "$OUT"
TIMES="$OUT/timings.tsv"; printf "step\tseconds\tstatus\n" > "$TIMES"
FORCE_FLAG=(); [[ "$KEEP" == "1" ]] || FORCE_FLAG=(--force)

step(){ # step <name> <cmd...>; records wall time and stops on failure
  local name="$1"; shift
  local t0 rc=0; t0=$(date +%s)
  echo "[e2e] ==== $name"
  "$@" > "$OUT/$name.log" 2>&1 || rc=$?
  printf "%s\t%s\t%s\n" "$name" "$(( $(date +%s) - t0 ))" "$([[ $rc -eq 0 ]] && echo ok || echo "fail($rc)")" >> "$TIMES"
  [[ $rc -eq 0 ]] || { echo "[e2e] $name failed (rc=$rc), log: $OUT/$name.log" >&2; tail -n 30 "$OUT/$name.log" >&2; exit "$rc"; }
}

need(){ [[ -s "$1" ]] || { echo "[e2e] expected output missing: $1" >&2; exit 4; }; }

# 1) simulate
step sim "$MAIN" run sim
REFDIR="$ROOT/refdata/testsets/$SET"
R1="$ROOT/SRA/reads/$SET/metagenome_R1.fastq.gz"; R2="$ROOT/SRA/reads/$SET/metagenome_R2.fastq.gz"
need "$R1"; need "$R2"; need "$REFDIR/refs/combined.fna"; need "$REFDIR/truth/provenance.tsv"

# 2) ingest via a one-row manifest
MANI="$OUT/manifest.tsv"
printf "sample\tfq1\tfq2\n%s\t%s\t%s\n" "$SET" "$(readlink -f "$R1")" "$(readlink -f "$R2")" > "$MANI"
step ingest   "$MAIN" run ingest --fastq "$MANI"
step qc       "$MAIN" run qc       --sample "$SET" "${FORCE_FLAG[@]}"
step assembly "$MAIN" run assembly --sample "$SET" "${FORCE_FLAG[@]}"
CONTIGS="$ROOT/SRA/assemblies/spades/$SET/contigs.fasta"; need "$CONTIGS"
if [[ "$SKIP_BINNING" == "0" ]]; then
  step binning "$MAIN" run binning --sample "$SET" "${FORCE_FLAG[@]}"
fi
step plasmids "$MAIN" run plasmids --sample "$SET" "${FORCE_FLAG[@]}"
UNION="$ROOT/SRA/plasmids/$SET/union/plasmids_raw.fasta"; need "$UNION"

# 3) score against the truth
EVAL_ENV="memamore_eval"
micromamba env list | awk '{print $1}' | grep -qx "$EVAL_ENV" \
  || micromamba create -y -q -n "$EVAL_ENV" -c conda-forge -c bioconda "minimap2=2.28" "python=3.11"
mm2(){ micromamba run -n "$EVAL_ENV" minimap2 -x asm10 -t "$CPUS" --secondary=no "$REFDIR/refs/combined.fna" "$1" 2>/dev/null; }
mm2 "$CONTIGS" > "$OUT/contigs_vs_refs.paf"
mm2 "$UNION"   > "$OUT/plasmids_vs_refs.paf"
micromamba run -n "$EVAL_ENV" python "$ROOT/tests/e2e/evaluate.py" \
  --refs "$REFDIR/refs/combined.fna" --provenance "$REFDIR/truth/provenance.tsv" \
  --contigs-paf "$OUT/contigs_vs_refs.paf" --pred-paf "$OUT/plasmids_vs_refs.paf" \
  --out-tsv "$OUT/metrics.tsv" --out-json "$OUT/metrics.json" > /dev/null

# MAG quality, if binning ran
Q="$ROOT/SRA/binning/checkm2/$SET/quality_report.tsv"
if [[ -s "$Q" ]]; then
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="Completeness")c=i; if($i=="Contamination")k=i}; next}
              {n++; if($c>=90&&$k<5)hq++; else if($c>=50&&$k<10)mq++}
              END{printf "bins\t%d\nhq_mimag\t%d\nmq_mimag\t%d\n", n, hq, mq}' "$Q" > "$OUT/mags.tsv"
fi

# 4) provenance for the paper: git commit + tool versions of every env
{
  printf "item\tvalue\n"
  printf "git_commit\t%s\n" "$(git -C "$ROOT" rev-parse HEAD 2>/dev/null || echo unknown)"
  printf "config\t%s\n" "$CFG"
  printf "date_utc\t%s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$OUT/versions.tsv"
for env in $(micromamba env list | awk 'NR>2 && $1!~/^#/ {print $1}'); do
  micromamba env export -n "$env" --explicit > "$OUT/env_${env}.lock" 2>/dev/null || true
done

echo "[e2e] done -> $OUT"
column -t -s $'\t' "$OUT/metrics.tsv"
[[ -s "$OUT/mags.tsv" ]] && column -t -s $'\t' "$OUT/mags.tsv"
column -t -s $'\t' "$TIMES"
