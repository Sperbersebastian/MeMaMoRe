#!/usr/bin/env python3
"""Stub HOTSPOT preprocessing: remember the input FASTA for hotspot.py."""
import argparse, os
p = argparse.ArgumentParser()
p.add_argument("--fasta"); p.add_argument("--database"); p.add_argument("--model_path")
a = p.parse_args()
os.makedirs("results", exist_ok=True)
open("results/.input", "w").write(a.fasta)
