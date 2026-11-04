#!/usr/bin/env python3
"""Build simulated 2x150 bp test samples from public NCBI references.

Profiles
  mini  (default) phage phiX174, phage lambda, E. coli plasmid F; ~150 kb, ~40x,
        ~2 MB of reads, one sample "mini". Exercises ingest/qc/assembly/geNomad.
        Expected: metaSPAdes recovers 3 contigs (~98, ~48.5, ~5.4 kb); geNomad
        calls F a plasmid and lambda/phiX viruses.
  mock  small mock community for testing every module, two samples "mockA" and
        "mockB" with different abundances (multi-sample summary, binning):
          S. aureus NCTC 8325 chromosome (2.8 Mb)   MAG, virulence genes, mecA-free
          H. pylori 26695 (1.7 Mb)                  MAG, cagA/vacA virulence genes
          M. pneumoniae M129 (0.8 Mb)               MAG, low coverage
          E. coli plasmid pEK499 (118 kb)           plasmid with CTX-M-15, TEM-1, OXA-1 ...
          phages T4, lambda, phiX174                viruses
        ~5.6 Mb, ~450k read pairs (~70 MB) per sample.

Usage: scripts/make_mini_testdata.py [OUTDIR] [--profile mini|mock]
       (default OUTDIR: test/<profile>)
Then:  bin/main.sh run ingest --fastq OUTDIR/manifest.tsv
"""
import argparse
import gzip
import os
import random
import urllib.request

URL = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
       "?db=nuccore&id={}&rettype=fasta&retmode=text")
READ_LEN, INSERT, ERR_RATE = 150, 350, 0.002

# profile -> {sample: {accession: coverage}}
PROFILES = {
    "mini": {
        "mini": {"NC_001422.1": 40, "NC_001416.1": 40, "NC_002483.1": 40},
    },
    "mock": {
        "mockA": {"NC_007795.1": 30, "NC_000915.1": 18, "NC_000912.1": 12,
                  "NC_013122.1": 45, "NC_000866.4": 30, "NC_001416.1": 35, "NC_001422.1": 50},
        "mockB": {"NC_007795.1": 12, "NC_000915.1": 30, "NC_000912.1": 20,
                  "NC_013122.1": 25, "NC_000866.4": 10, "NC_001416.1": 60, "NC_001422.1": 20},
    },
}


def fetch(accs, refs):
    if not os.path.isfile(refs):
        with urllib.request.urlopen(URL.format(",".join(accs))) as r, open(refs + ".part", "wb") as f:
            f.write(r.read())
        os.replace(refs + ".part", refs)
    seqs, name = {}, None
    for line in open(refs):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif line:
            seqs[name].append(line.upper())
    seqs = {k: "".join(v) for k, v in seqs.items()}
    missing = set(accs) - set(seqs)
    if missing:
        raise SystemExit(f"[testdata] missing in {refs}: {sorted(missing)} (delete it to re-download)")
    return seqs


def simulate(sample, cov, seqs, out, seed):
    rnd = random.Random(seed)
    comp = str.maketrans("ACGTN", "TGCAN")

    def rc(s):
        return s.translate(comp)[::-1]

    def mutate(s):
        # draw error positions instead of testing every base (fast for big samples)
        s = list(s)
        pos = 0
        while True:
            pos += int(rnd.expovariate(ERR_RATE)) + 1
            if pos > len(s):
                return "".join(s)
            s[pos - 1] = rnd.choice("ACGT")

    # Varied Phred+33 qualities (Q20-Q40). A constant "I" string is ambiguous
    # between Phred+33 and Phred+64, and SPAdes' error correction then aborts
    # with "Failed to determine offset".
    quals = ["".join(chr(33 + min(40, max(20, int(rnd.gauss(36 - 6 * i / READ_LEN, 3)))))
                     for i in range(READ_LEN)) for _ in range(64)]
    r1 = os.path.join(out, f"{sample}_R1.fastq.gz")
    r2 = os.path.join(out, f"{sample}_R2.fastq.gz")
    n = 0
    with gzip.open(r1, "wt", compresslevel=3) as f1, gzip.open(r2, "wt", compresslevel=3) as f2:
        for acc, c in cov.items():
            seq = seqs[acc]
            for _ in range(len(seq) * c // (2 * READ_LEN)):
                ins = max(2 * READ_LEN, int(rnd.gauss(INSERT, 30)))
                p = rnd.randint(0, len(seq) - ins)
                frag = seq[p:p + ins]
                if rnd.random() < 0.5:
                    frag = rc(frag)
                n += 1
                f1.write(f"@r{n}/1\n{mutate(frag[:READ_LEN])}\n+\n{rnd.choice(quals)}\n")
                f2.write(f"@r{n}/2\n{mutate(rc(frag)[:READ_LEN])}\n+\n{rnd.choice(quals)}\n")
    print(f"[testdata] {sample}: {n} read pairs from "
          + ", ".join(f"{a} ({len(seqs[a])} bp, {c}x)" for a, c in cov.items()))
    return r1, r2


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("outdir", nargs="?")
    ap.add_argument("--profile", choices=sorted(PROFILES), default="mini")
    a = ap.parse_args()
    out = os.path.abspath(a.outdir or f"test/{a.profile}")
    os.makedirs(out, exist_ok=True)

    samples = PROFILES[a.profile]
    accs = sorted({acc for cov in samples.values() for acc in cov})
    seqs = fetch(accs, os.path.join(out, "refs.fasta"))

    rows = []
    for i, (sample, cov) in enumerate(samples.items(), 1):
        r1, r2 = simulate(sample, cov, seqs, out, seed=i)
        rows.append(f"{sample}\t{r1}\t{r2}")
    with open(os.path.join(out, "manifest.tsv"), "w") as m:
        m.write("sample\tfq1\tfq2\n" + "\n".join(rows) + "\n")
    print(f"[testdata] manifest: {os.path.join(out, 'manifest.tsv')}")


if __name__ == "__main__":
    main()
