#!/usr/bin/env python3
"""Stub HOTSPOT preprocessing: remember the input FASTA for hotspot.py."""
import argparse
import os

p = argparse.ArgumentParser()
for opt in ("--fasta", "--database", "--model_path", "--threads", "--len"):
    p.add_argument(opt)
p.add_argument("--midfolder", default="temp")
a = p.parse_args()
os.makedirs(a.midfolder, exist_ok=True)
with open(os.path.join(a.midfolder, "input.txt"), "w") as f:
    f.write(a.fasta)
