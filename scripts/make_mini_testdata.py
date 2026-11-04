#!/usr/bin/env python3
"""Build the tiny 'mini' test sample: 2x150 bp reads simulated from phage phiX174,
phage lambda and E. coli plasmid F (~150 kb total, ~40x, ~2 MB of reads).

Expected results: metaSPAdes recovers 3 contigs (~98, ~48.5, ~5.4 kb);
geNomad calls F a plasmid and lambda/phiX viruses.

Usage: scripts/make_mini_testdata.py [OUTDIR]   (default: test/mini)
Then:  bin/main.sh run ingest --fastq OUTDIR/manifest.tsv
"""
import gzip, os, random, sys, urllib.request

ACCS = ["NC_001422.1", "NC_001416.1", "NC_002483.1"]  # phiX174, lambda, plasmid F
URL = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
       "?db=nuccore&id={}&rettype=fasta&retmode=text")
READ_LEN, INSERT, COVERAGE, ERR_RATE, SEED = 150, 350, 40, 0.002, 1

def main():
    out = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "test/mini")
    os.makedirs(out, exist_ok=True)
    refs = os.path.join(out, "refs.fasta")
    if not os.path.isfile(refs):
        with urllib.request.urlopen(URL.format(",".join(ACCS))) as r, open(refs, "wb") as f:
            f.write(r.read())

    seqs, name = {}, None
    for line in open(refs):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        elif line:
            seqs[name].append(line.upper())
    seqs = {k: "".join(v) for k, v in seqs.items()}

    rnd = random.Random(SEED)
    comp = str.maketrans("ACGTN", "TGCAN")
    rc = lambda s: s.translate(comp)[::-1]
    def mutate(s):
        return "".join(rnd.choice("ACGT") if rnd.random() < ERR_RATE else c for c in s)

    r1 = os.path.join(out, "mini_R1.fastq.gz"); r2 = os.path.join(out, "mini_R2.fastq.gz")
    n, qual = 0, "I" * READ_LEN
    with gzip.open(r1, "wt") as f1, gzip.open(r2, "wt") as f2:
        for seq in seqs.values():
            for _ in range(len(seq) * COVERAGE // (2 * READ_LEN)):
                ins = max(2 * READ_LEN, int(rnd.gauss(INSERT, 30)))
                p = rnd.randint(0, len(seq) - ins)
                frag = seq[p:p + ins]
                if rnd.random() < 0.5:
                    frag = rc(frag)
                n += 1
                f1.write(f"@r{n}/1\n{mutate(frag[:READ_LEN])}\n+\n{qual}\n")
                f2.write(f"@r{n}/2\n{mutate(rc(frag)[:READ_LEN])}\n+\n{qual}\n")

    with open(os.path.join(out, "manifest.tsv"), "w") as m:
        m.write(f"sample\tfq1\tfq2\nmini\t{r1}\t{r2}\n")
    print(f"[mini] {n} read pairs from {', '.join(f'{k} ({len(v)} bp)' for k, v in seqs.items())}")
    print(f"[mini] manifest: {os.path.join(out, 'manifest.tsv')}")

if __name__ == "__main__":
    main()
