#!/usr/bin/env python3
"""Binning groups of all samples.

Usage: sample_groups.py ROOT            -> TSV "sample<TAB>group" (group may be empty)
       sample_groups.py ROOT SAMPLE     -> members of SAMPLE's group, one per line
                                           (just SAMPLE when it has no group)

Groups come from the 'group' column of SRA/reads/manifest.tsv (set at ingest);
a 'group' column in config/samples.tsv (key column 'sample' or 'sample_id')
overrides it. Samples of one group are binned with each other's coverage;
assembly always stays per sample. Standard library only.
"""
import csv
import os
import sys


def read_rows(path):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def groups(root):
    out = {}
    for r in read_rows(os.path.join(root, "SRA", "reads", "manifest.tsv")):
        if r.get("sample"):
            out[r["sample"]] = (r.get("group") or "").strip()
    for r in read_rows(os.path.join(root, "config", "samples.tsv")):
        s = (r.get("sample") or r.get("sample_id") or "").strip()
        g = (r.get("group") or "").strip()
        if s and g and s in out:
            out[s] = g
    return out


def main():
    root = sys.argv[1]
    g = groups(root)
    if len(sys.argv) > 2:
        s = sys.argv[2]
        grp = g.get(s, "")
        members = sorted(x for x, y in g.items() if grp and y == grp) or [s]
        print("\n".join(members))
        return
    for s in sorted(g):
        print(f"{s}\t{g[s]}")


if __name__ == "__main__":
    main()
