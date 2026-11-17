"""Unit tests for scripts/plasmid_hosts.py (plasmid host by abundance correlation)."""
import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "plasmid_hosts.py"
SAMPLES = [f"s{i}" for i in range(1, 7)]


def write_tsv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def test_host_follows_coverage_not_tpm(tmp_path):
    """A plasmid that follows host A is assigned to A although a large chromosome
    fragment of B in the plasmid pool makes its TPM profile follow B inversely."""
    sra = tmp_path / "SRA"
    host_a = [10, 20, 30, 40, 50, 60]           # % of reads
    host_b = [5, 80, 10, 90, 20, 95]            # independent of A
    nreads = [1e6, 3e6, 2e6, 1.5e6, 4e6, 2.5e6]  # sequencing depth differs per sample
    gfa = sra / "plasmids" / "_global" / "cluster" / "plasmids_derep.global.fasta"
    gfa.parent.mkdir(parents=True)
    gfa.write_text(">s1|genomad|plasmidA\n" + "A" * 5000 + "\n>s1|genomad|fragB\n" + "C" * 200000 + "\n")
    for s, a, b, n in zip(SAMPLES, host_a, host_b, nreads):
        hj = sra / "qc" / "hostile" / s / "hostile.json"
        hj.parent.mkdir(parents=True)
        hj.write_text(json.dumps([{"reads_in": n, "reads_out": n}]))
        # mean coverage scales with abundance x depth; TPM is relative within the pool
        cov_p, cov_f = a * 2 * n / 1e6, b * n / 1e6
        rpk_p, rpk_f = cov_p, cov_f
        tpm_p, tpm_f = 1e6 * rpk_p / (rpk_p + rpk_f), 1e6 * rpk_f / (rpk_p + rpk_f)
        reads = f"plasmids_derep.global.fasta/{s}_clean_1.fastq.gz"
        write_tsv(tmp_path / "catalog" / "plasmid_abundance" / f"{s}.tsv",
                  ["Contig", f"{reads} TPM", f"{reads} Mean", f"{reads} Covered Fraction"],
                  [["s1|genomad|plasmidA", tpm_p, cov_p, 1], ["s1|genomad|fragB", tpm_f, cov_f, 1]])
    hdr = ["catalog_genome", "gtdb_classification"] + [f"{s}_rel_abundance_pct" for s in SAMPLES]
    write_tsv(sra / "summary" / "mag_catalog.tsv", hdr,
              [["magA", "s__A", *host_a], ["magB", "s__B", *host_b]])
    subprocess.run([sys.executable, str(SCRIPT), str(tmp_path), str(tmp_path / "catalog"),
                    "6", "0.7", "0.05", *SAMPLES], check=True, capture_output=True)
    rows = {r["plasmid"]: r for r in csv.DictReader(open(sra / "summary" / "plasmid_hosts.tsv"), delimiter="\t")}
    assert rows["s1|genomad|plasmidA"]["corr_host_genome"] == "magA"
    assert rows["s1|genomad|fragB"]["corr_host_genome"] == "magB"
