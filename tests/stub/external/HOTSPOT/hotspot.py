#!/usr/bin/env python3
"""Stub HOTSPOT prediction: one host lineage row per input contig."""
fa = open("results/.input").read().strip()
ids = [l[1:].split()[0] for l in open(fa) if l.startswith(">")]
with open("results/host_lineage.tsv", "w") as f:
    f.write("contig\tphylum\tclass\torder\tfamily\tgenus\n")
    for i in ids:
        f.write(f"{i}\tPseudomonadota\tGammaproteobacteria\tEnterobacterales\tEnterobacteriaceae\tEscherichia\n")
