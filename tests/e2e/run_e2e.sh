#!/usr/bin/env bash
# End-to-end test on a simulated community with known composition.
#
#   tests/e2e/run_e2e.sh [--config config/sim.yaml] [--cpus N] [--skip-binning] [--keep]
#   tests/e2e/run_e2e.sh --external DIR [--name NAME] ...   # reads and truth from elsewhere
#
# --external DIR skips the simulation and uses DIR/manifest.tsv (sample, fq1, fq2[, group]),
# DIR/refs/combined.fna and DIR/truth/provenance.tsv (type, accession, ...), e.g. a published
# benchmark data set prepared by MeMaMoRe-paper/benchmark/external/*/prepare.sh.
#
# Steps: simulate reads (modules/sim; one or more samples, all in one group) -> ingest ->
#        per sample qc -> assembly -> [binning] -> plasmids -> align contigs and plasmid
#        predictions to the known references -> metrics (per tool, union, >= 2 and >= 3 tools),
#        MetaQUAST, MAG quality.
# Results: tests/e2e/results/<SET>/{metrics.tsv,timings.tsv,versions.tsv,provenance.tsv,weights.tsv}
#          plus one folder per sample (metrics.json, metaquast/, mags.tsv, PAF files)
#
# Needs micromamba plus the databases the modules use (geNomad, PLASMe, CheckM2, GTDB-Tk
# for binning). Not run in CI: it downloads references and takes hours on a workstation.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CFG="$ROOT/config/sim.yaml"
CPUS="${CPUS:-8}"
SKIP_BINNING=0
KEEP=0
EXTERNAL=""
NAME=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CFG="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"; shift 2 ;;
    --cpus) CPUS="$2"; shift 2 ;;
    --skip-binning) SKIP_BINNING=1; shift ;;
    --keep) KEEP=1; shift ;;
    --external) EXTERNAL="$(cd "$2" && pwd)"; shift 2 ;;
    --name) NAME="$2"; shift 2 ;;
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
if [[ -n "$EXTERNAL" ]]; then
  SET="${NAME:-$(basename "$EXTERNAL")}"
else
  SET="$(grep -E '^set_name:' "$CFG" | head -n1 | sed -E 's/^set_name:[[:space:]]*"?([^"#]*)"?.*/\1/' | xargs)"
fi
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

# 1) simulate (one or more samples, all in group $SET; see scripts/sim_plan.py)
if [[ -n "$EXTERNAL" ]]; then
  REFDIR="$EXTERNAL"; SIMMANI="$EXTERNAL/manifest.tsv"
  [[ -s "$REFDIR/truth/weights.tsv" ]] || printf "sample\taccession\tweight\n" > "$OUT/weights.tsv"
else
  step sim "$MAIN" run sim
  REFDIR="$ROOT/refdata/testsets/$SET"; SIMMANI="$REFDIR/sim/manifest.tsv"
fi
need "$SIMMANI"; need "$REFDIR/refs/combined.fna"; need "$REFDIR/truth/provenance.tsv"
mapfile -t SAMPLES < <(tail -n +2 "$SIMMANI" | cut -f1)
[[ ${#SAMPLES[@]} -gt 0 ]] || { echo "[e2e] no samples in $SIMMANI" >&2; exit 4; }
cp "$REFDIR/truth/provenance.tsv" "$OUT/"
[[ -s "$REFDIR/truth/weights.tsv" ]] && cp "$REFDIR/truth/weights.tsv" "$OUT/"

# 2) ingest the simulation manifest (sample, fq1, fq2, group), then every module per sample
step ingest "$MAIN" run ingest --fastq "$SIMMANI"
for S in "${SAMPLES[@]}"; do
  step "qc_$S"       "$MAIN" run qc       --sample "$S" "${FORCE_FLAG[@]}"
  step "assembly_$S" "$MAIN" run assembly --sample "$S" "${FORCE_FLAG[@]}"
  need "$ROOT/SRA/assemblies/spades/$S/contigs.fasta"
  if [[ "$SKIP_BINNING" == "0" ]]; then
    step "binning_$S" "$MAIN" run binning --sample "$S" "${FORCE_FLAG[@]}"
  fi
  step "plasmids_$S" "$MAIN" run plasmids --sample "$S" "${FORCE_FLAG[@]}"
  need "$ROOT/SRA/plasmids/$S/union/plasmids_raw.fasta"
done

# 3) score every sample against the truth
EVAL_ENV="memamore_eval"
micromamba env list | awk '{print $1}' | grep -qx "$EVAL_ENV" \
  || micromamba create -y -q -n "$EVAL_ENV" -c conda-forge -c bioconda "minimap2=2.28" "quast=5.2" "python=3.11"
mm2(){ micromamba run -n "$EVAL_ENV" minimap2 -x asm10 -t "$CPUS" --secondary=no "$1" "$2" 2>/dev/null; }
for S in "${SAMPLES[@]}"; do
  SO="$OUT/$S"; mkdir -p "$SO"
  CONTIGS="$ROOT/SRA/assemblies/spades/$S/contigs.fasta"
  UNION="$ROOT/SRA/plasmids/$S/union/plasmids_raw.fasta"
  mm2 "$REFDIR/refs/combined.fna" "$CONTIGS" > "$SO/contigs_vs_refs.paf"
  mm2 "$REFDIR/refs/combined.fna" "$UNION"   > "$SO/plasmids_vs_refs.paf"
  # all-vs-all of the predictions: which tools agree (consensus rows min2/min3)
  micromamba run -n "$EVAL_ENV" minimap2 -x asm10 -X -t "$CPUS" "$UNION" "$UNION" 2>/dev/null > "$SO/plasmids_self.paf"
  micromamba run -n "$EVAL_ENV" python "$ROOT/tests/e2e/evaluate.py" \
    --refs "$REFDIR/refs/combined.fna" --provenance "$REFDIR/truth/provenance.tsv" \
    --contigs-paf "$SO/contigs_vs_refs.paf" --pred-paf "$SO/plasmids_vs_refs.paf" \
    --self-paf "$SO/plasmids_self.paf" \
    --out-tsv "$SO/metrics.tsv" --out-json "$SO/metrics.json" > /dev/null
  # assembly quality against the references
  micromamba run -n "$EVAL_ENV" metaquast.py -o "$SO/metaquast" -r "$REFDIR/refs/combined.fna" \
    -t "$CPUS" --min-contig 1000 --silent "$CONTIGS" > "$SO/metaquast.log" 2>&1 \
    || echo "[e2e] metaquast failed for $S (see $SO/metaquast.log)"
  # MAG quality, if binning ran
  Q="$ROOT/SRA/binning/checkm2/$S/quality_report.tsv"
  if [[ -s "$Q" ]]; then
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="Completeness")c=i; if($i=="Contamination")k=i}; next}
                {n++; if($c>=90&&$k<5)hq++; else if($c>=50&&$k<10)mq++}
                END{printf "bins\t%d\nhq_mimag\t%d\nmq_mimag\t%d\n", n, hq, mq}' "$Q" > "$SO/mags.tsv"
  fi
done

# one table over all samples: sample + per-tool metrics
first=1; : > "$OUT/metrics.tsv"
for S in "${SAMPLES[@]}"; do
  if [[ $first -eq 1 ]]; then sed -n 1p "$OUT/$S/metrics.tsv" | sed 's/^/sample\t/' >> "$OUT/metrics.tsv"; first=0; fi
  tail -n +2 "$OUT/$S/metrics.tsv" | sed "s/^/$S\t/" >> "$OUT/metrics.tsv"
done

# 4) provenance for the paper: git commit + tool versions of every env
{
  printf "item\tvalue\n"
  printf "git_commit\t%s\n" "$(git -C "$ROOT" rev-parse HEAD 2>/dev/null || echo unknown)"
  printf "config\t%s\n" "${EXTERNAL:-$CFG}"
  printf "date_utc\t%s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$OUT/versions.tsv"
for env in $(micromamba env list | awk 'NR>2 && $1!~/^#/ {print $1}'); do
  micromamba env export -n "$env" --explicit > "$OUT/env_${env}.lock" 2>/dev/null || true
done

echo "[e2e] done -> $OUT"
show(){ if command -v column >/dev/null; then column -t -s $'\t' "$1"; else cat "$1"; fi; }
show "$OUT/metrics.tsv"
show "$TIMES"
