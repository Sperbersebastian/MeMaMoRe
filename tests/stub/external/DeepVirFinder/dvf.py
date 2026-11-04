#!/usr/bin/env python3
"""Stub DeepVirFinder: dvf.py -i IN -o OUT -l LEN -c N -> OUT/<IN>_gt<LEN>bp_dvfpred.txt
(name len score pvalue); contigs of the toy phages (NODE_2, NODE_6) score as viral."""
import argparse
import os

p = argparse.ArgumentParser()
p.add_argument("-i"); p.add_argument("-o"); p.add_argument("-l", type=int, default=1); p.add_argument("-c")
a = p.parse_args()
recs, cur = [], None
for ln in open(a.i):
    if ln.startswith(">"):
        cur = [ln[1:].split()[0], 0]
        recs.append(cur)
    elif cur:
        cur[1] += len(ln.strip())
os.makedirs(a.o, exist_ok=True)
with open(os.path.join(a.o, f"{os.path.basename(a.i)}_gt{a.l}bp_dvfpred.txt"), "w") as f:
    f.write("name\tlen\tscore\tpvalue\n")
    for name, n in recs:
        if n >= a.l:
            viral = name.startswith(("NODE_2_", "NODE_6_"))
            f.write(f"{name}\t{n}\t{0.99 if viral else 0.05}\t{0.001 if viral else 0.8}\n")
