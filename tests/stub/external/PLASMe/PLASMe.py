#!/usr/bin/env python3
"""Stub PLASMe: PLASMe.py <in.fasta> <out.fasta> [-d DB -m MODE -t N] -> plasmid contigs + report."""
import sys
inp, out = sys.argv[1], sys.argv[2]
recs, cur = [], None
for line in open(inp):
    if line.startswith(">"):
        cur = [line[1:].split()[0], []]
        recs.append(cur)
    elif cur:
        cur[1].append(line.strip())
keep = [(c, "".join(s)) for c, s in recs if c.startswith("NODE_1_")]
with open(out, "w") as f:
    for c, s in keep:
        f.write(f">{c}\n{s}\n")
with open(out.rsplit(".", 1)[0] + "_report.tsv", "w") as f:
    f.write("contig\tlength\treference\torder\tevidence\tscore\tamb_region\n")
    for c, s in keep:
        f.write(f"{c}\t{len(s)}\tstub_ref\tEnterobacterales\tBLASTn\t0.98\t-\n")
