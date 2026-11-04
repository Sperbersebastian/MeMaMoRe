#!/usr/bin/env bash
# Stub run: the whole pipeline with fake tools, no envs, no databases.
#
#   tests/stub/run.sh [module ...]     (default: ingest qc assembly binning plasmids viruses args)
#
# Runs in a throwaway ROOT (STUB_WORK, default: a new temp dir) that symlinks
# the code from this repo, so real SRA/, refdata/ and logs/ are never touched.
# A fake `micromamba` (tests/stub/bin) pretends every env exists and runs
# commands directly; fake tools next to it write correctly named outputs.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
STUB="$REPO/tests/stub"
WORK="${STUB_WORK:-$(mktemp -d -t memamo_stub.XXXXXX)}"
SAMPLE="${STUB_SAMPLE:-stub}"
MODULES=("$@"); [[ ${#MODULES[@]} -gt 0 ]] || MODULES=(ingest qc assembly binning plasmids viruses args)
# extra config for every module, e.g. STUB_SET="modules.qc.host_removal=human modules.plasmids.consensus=min2"
SETS=(--set global.threads="${CPUS:-2}")
for kv in ${STUB_SET:-}; do SETS+=(--set "$kv"); done

mkdir -p "$WORK"
for d in bin modules envs config scripts gui; do ln -sfn "$REPO/$d" "$WORK/$d"; done

export STUB_PYTHON="${STUB_PYTHON:-$(command -v python3)}"
export MAMBA_ROOT_PREFIX="$STUB"          # main.sh/env.sh prepend $MAMBA_ROOT_PREFIX/bin
export PATH="$STUB/bin:/usr/local/bin:/usr/bin:/bin"
export MEMAMO_STUB=1 CPUS="${CPUS:-2}" GTDBTK_FORCE=1 PLASMER_FORCE=1   # stub: exercise RAM-gated steps
# comma-decimal locale (if installed) so CI catches awk/sort number parsing bugs
locale -a 2>/dev/null | grep -qi '^de_DE.utf' && export LC_ALL='' LANG=de_DE.UTF-8 LC_NUMERIC=de_DE.UTF-8
unset ROOT GTDBTK_DATA_PATH REFDATA_BASE VIBRANT_DB

# --- placeholder external tools / reference data (only what scripts test for) ---
fake(){ mkdir -p "$(dirname "$WORK/$1")"; [[ -s "$WORK/$1" ]] || echo "stub" > "$WORK/$1"; }
fake external/binny/Snakefile
fake external/binny/config/binny_mantis.cfg
printf 'NAME  stub_tf\nACC   TIGR00001\nLENG  100\n' > "$WORK/stub.hmm"
for k in tf pf; do
  mkdir -p "$WORK/external/binny/database/hmms/checkm_$k"
  cp "$WORK/stub.hmm" "$WORK/external/binny/database/hmms/checkm_$k/checkm_filtered_$k.hmm"
done
fake external/MAGScoT/MAGScoT.R
fake external/MAGScoT/hmm/gtdbtk_rel207_tigrfam.hmm
fake external/MAGScoT/hmm/gtdbtk_rel207_Pfam-A.hmm
mkdir -p "$WORK/refdata/gtdbtk/release226/markers" "$WORK/refdata/mobsuite_db" "$WORK/refdata/plasme/DB"
fake refdata/viralverify/pfam/Pfam-A.hmm
fake refdata/genomad/genomad_db/version.txt
mkdir -p "$WORK/refdata/iphop/Aug_2023_pub_rw"
fake refdata/gunc/gunc_db_progenomes2.1.dmnd
fake refdata/plasmer/.extracted
fake refdata/mobileog/mobileOGdb_90-2.0.fasta
fake refdata/virsorter2/db/Done_all_setup
for h in Pfam-A_v32 KEGG_profiles_prokaryotes VOGDB94_phage; do for x in h3f h3i h3m h3p; do fake "refdata/vibrant/db/databases/$h.HMM.$x"; done; done
fake refdata/vibrant/db/files/VIBRANT_names.tsv
fake refdata/checkv/checkv-db-v1.5/genome_db/checkv_reps.faa
fake refdata/checkv/checkv-db-v1.5/genome_db/checkv_reps.dmnd
mkdir -p "$WORK/refdata/amrplusplus/megares_v3" "$WORK/refdata/card/localDB"
printf '>MEG_1|Drugs|betalactams|Class_A_betalactamases|TEM\nACGTACGTACGT\n>MEG_2|Drugs|Tetracyclines|Tetracycline_efflux|TETA\nACGTACGTACGT\n' \
  > "$WORK/refdata/amrplusplus/megares_v3/megares_database_v3.00.fasta"
printf 'header,class,mechanism,group\nMEG_1|Drugs|betalactams|Class_A_betalactamases|TEM,betalactams,Class A betalactamases,TEM\n' \
  > "$WORK/refdata/amrplusplus/megares_v3/megares_annotations_v3.00.csv"
# small placeholder scripts for tools that are run from external/ checkouts
mkdir -p "$WORK/external"
cp -rn "$STUB/external/." "$WORK/external/"

# --- tiny paired-end reads for two samples in one binning group ---------------
# (all samples of a test run share a group, so multi-sample binning is always tested)
SAMPLES=("$SAMPLE" "${SAMPLE}2")
mkdir -p "$WORK/input"
"$STUB_PYTHON" - "$WORK/input" "${SAMPLES[@]}" <<'PY'
import gzip, random, sys
out = sys.argv[1]
for k, s in enumerate(sys.argv[2:]):
    rnd = random.Random(1 + k)
    for mate in (1, 2):
        with gzip.open(f"{out}/{s}_{mate}.fastq.gz", "wt") as f:
            for i in range(200):
                seq = "".join(rnd.choice("ACGT") for _ in range(150))
                f.write(f"@r{i}/{mate}\n{seq}\n+\n{'I' * 150}\n")
PY
{
  printf 'sample\tfq1\tfq2\tgroup\n'
  for s in "${SAMPLES[@]}"; do
    printf '%s\t%s\t%s\t%s\n' "$s" "$WORK/input/${s}_1.fastq.gz" "$WORK/input/${s}_2.fastq.gz" "stubgroup"
  done
} > "$WORK/input/manifest.tsv"

echo "[stub] ROOT=$WORK samples=${SAMPLES[*]} (group stubgroup) modules=${MODULES[*]}"
for m in "${MODULES[@]}"; do
  echo "[stub] ===== $m ====="
  case "$m" in
    ingest)  "$WORK/bin/main.sh" run ingest --fastq "$WORK/input/manifest.tsv" ;;
    binning) "$WORK/bin/main.sh" run binning --sample "$SAMPLE" "${SETS[@]}" ;;   # bins the whole group
    *)       for s in "${SAMPLES[@]}"; do "$WORK/bin/main.sh" run "$m" --sample "$s" "${SETS[@]}"; done ;;
  esac
done
echo "[stub] ===== catalog ====="
if [[ " ${MODULES[*]} " == *" binning "* ]]; then "$WORK/bin/main.sh" run catalog "${SETS[@]}"; fi
echo "[stub] ===== summary ====="
"$WORK/bin/main.sh" summary
echo "[stub] all modules finished: ${MODULES[*]}"
echo "[stub] outputs under $WORK/SRA"
