#!/usr/bin/env python3
"""Stub HOTSPOT prediction: one host lineage row per input contig."""
import argparse
import os

p = argparse.ArgumentParser()
for opt in ("--model_path", "--threads", "--device", "--batch_size"):
    p.add_argument(opt)
p.add_argument("--midfolder", default="temp")
p.add_argument("--out", default="results")
a = p.parse_args()
fa = open(os.path.join(a.midfolder, "input.txt")).read().strip()
ids = [ln[1:].split()[0] for ln in open(fa) if ln.startswith(">")]
os.makedirs(a.out, exist_ok=True)
with open(os.path.join(a.out, "host_lineage.tsv"), "w") as f:
    f.write("Contig\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n")   # real HOTSPOT header
    for i in ids:
        f.write(f"{i}\tPseudomonadota\tGammaproteobacteria\tEnterobacterales\tEnterobacteriaceae\tEscherichia\t-\n")
