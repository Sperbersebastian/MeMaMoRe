#!/usr/bin/env python3
"""Fake download and read-simulation tools for tests/sim/run.sh (no network, no ART).

efetch           -> one random sequence named after -id (length from the accession's last digit)
datasets         -> a zip holding a genome package with a chromosome and one plasmid record
art_illumina     -> paired FASTQ files <prefix>1.fq / <prefix>2.fq with reads proportional to -f
"""
import os
import random
import sys
import zipfile

TOOL = os.path.basename(sys.argv[0])
ARGS = sys.argv[1:]


def opt(flag, default=None):
    return ARGS[ARGS.index(flag) + 1] if flag in ARGS else default


def seq(n, rnd):
    return "".join(rnd.choice("ACGT") for _ in range(n))


if TOOL == "efetch":
    acc = opt("-id")
    rnd = random.Random(acc)
    n = 2000 + 500 * (sum(map(ord, acc)) % 7)
    sys.stdout.write(f">{acc} stub sequence\n{seq(n, rnd)}\n")
elif TOOL == "datasets":
    acc, out = ARGS[ARGS.index("accession") + 1], opt("--filename")
    rnd = random.Random(acc)
    with zipfile.ZipFile(out, "w") as z:
        z.writestr(f"ncbi_dataset/data/{acc}/{acc}_genomic.fna",
                   f">NZ_CHR{acc[4:10]}.1 Stub bacterium chromosome, complete genome\n{seq(6000, rnd)}\n"
                   f">NZ_PL{acc[4:10]}.1 Stub bacterium plasmid pStub, complete sequence\n{seq(2500, rnd)}\n")
elif TOOL == "art_illumina":
    if "--help" in ARGS:
        sys.exit(0)
    fold, ln, prefix = float(opt("-f")), int(opt("-l")), opt("-o")
    rnd = random.Random(opt("-rs"))
    n = max(1, int(fold * 10))
    for mate in (1, 2):
        with open(f"{prefix}{mate}.fq", "w") as fh:
            for i in range(n):
                fh.write(f"@{os.path.basename(prefix)}{i}/{mate}\n{seq(ln, rnd)}\n+\n{'I' * ln}\n")
else:
    sys.exit(f"unknown fake tool {TOOL}")
