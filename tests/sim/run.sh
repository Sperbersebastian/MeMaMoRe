#!/usr/bin/env bash
# Simulation module test with fake download/ART tools (no network, no envs).
#
#   tests/sim/run.sh
#
# Runs bin/main.sh run sim twice in a throwaway ROOT: the legacy single-sample config
# (config/sim.yaml) and an extended config with a genome assembly, host-linked plasmids
# and generated per-sample abundance profiles. Checks read files, manifest, truth table
# and that a second run gives identical reads.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
WORK="${SIM_WORK:-$(mktemp -d -t memamo_sim.XXXXXX)}"
mkdir -p "$WORK"
for d in bin modules envs config scripts; do ln -sfn "$REPO/$d" "$WORK/$d"; done

export STUB_PYTHON="${STUB_PYTHON:-$(command -v python3)}"
export MAMBA_ROOT_PREFIX="$REPO/tests/stub"
export PATH="$REPO/tests/sim/bin:$REPO/tests/stub/bin:/usr/local/bin:/usr/bin:/bin"
export MEMAMO_STUB=1
unset ROOT

fail(){ echo "[sim-test] FAIL: $*" >&2; exit 1; }

# 1) legacy config: one sample, historical paths
CFG="$REPO/config/sim.yaml" "$WORK/bin/main.sh" run sim > "$WORK/legacy.log" 2>&1 || { cat "$WORK/legacy.log"; fail "legacy run"; }
[[ -s "$WORK/SRA/reads/miniset/metagenome_R1.fastq.gz" ]] || fail "legacy reads missing"
[[ "$(tail -n +2 "$WORK/refdata/testsets/miniset/sim/manifest.tsv" | wc -l)" -eq 1 ]] || fail "legacy manifest should have 1 sample"
grep -q $'^host\tNC_000913.3\tminiset\t' "$WORK/refdata/testsets/miniset/truth/provenance.tsv" || fail "legacy truth"

# 2) extended config: assembly with its own plasmid, host-linked plasmid, 3 generated samples
cat > "$WORK/ext.yaml" <<'YAML'
set_name: "simtest"
read_len: 100
n_reads: 2000
seed: 7
references:
  - {acc: GCF_000005845.2, type: host, abundance: 0.5}
  - {acc: NC_002516.2, type: host, abundance: 0.3}
  - {acc: NC_013353.1, type: plasmid, host: NC_002516.2, abundance: 0.1}
  - {acc: NC_001416.1, type: phage, abundance: 0.1}
profiles: {n: 3, sigma: 1.0, copy_sigma: 0.1}
YAML
run_ext(){ CFG="$WORK/ext.yaml" "$WORK/bin/main.sh" run sim > "$WORK/ext.log" 2>&1 || { cat "$WORK/ext.log"; fail "extended run"; }; }
run_ext
T="$WORK/refdata/testsets/simtest"
[[ "$(tail -n +2 "$T/sim/manifest.tsv" | wc -l)" -eq 3 ]] || fail "expected 3 samples in manifest"
awk -F'\t' 'NR>1 && $4!="simtest"{exit 1}' "$T/sim/manifest.tsv" || fail "all test samples must share group simtest"
for s in simtest_s1 simtest_s2 simtest_s3; do
  [[ -s "$WORK/SRA/reads/$s/metagenome_R2.fastq.gz" ]] || fail "reads missing for $s"
done
grep -q $'^plasmid\tNZ_PL000005.1\tsimtest\tGCF_000005845.2\t' "$T/truth/provenance.tsv" || fail "assembly plasmid not typed plasmid with host"
grep -q $'^host\tNZ_CHR000005.1\tsimtest\t\t' "$T/truth/provenance.tsv" || fail "assembly chromosome not typed host"
grep -q $'^plasmid\tNC_013353.1\tsimtest\tNC_002516.2\t' "$T/truth/provenance.tsv" || fail "host link of listed plasmid lost"
grep -c '^>' "$T/refs/combined.fna" | grep -qx 5 || fail "combined.fna should hold 5 records"
# a host-linked plasmid follows its host across samples (ratio varies only by copy_sigma)
python3 - "$T/truth/weights.tsv" <<'PY' || fail "plasmid/host weights do not co-vary"
import csv, sys, statistics as st
w = {}
for r in csv.DictReader(open(sys.argv[1]), delimiter="\t"):
    w.setdefault(r["sample"], {})[r["accession"]] = float(r["weight"])
ratio = [s["NC_013353.1"] / s["NC_002516.2"] for s in w.values()]
host = [s["NC_002516.2"] for s in w.values()]
assert st.pstdev(ratio) / st.mean(ratio) < st.pstdev(host) / st.mean(host), (ratio, host)
PY
# determinism: same config -> identical reads
sum1="$(zcat "$T"/sim/simtest_s2/metagenome_R1.fastq.gz | md5sum)"
run_ext
sum2="$(zcat "$T"/sim/simtest_s2/metagenome_R1.fastq.gz | md5sum)"
[[ "$sum1" == "$sum2" ]] || fail "reads differ between identical runs"

echo "[sim-test] ok ($WORK)"
