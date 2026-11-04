#!/usr/bin/env python3
"""Plasmid consensus filter for union_cluster.

Usage: plasmid_consensus.py IN.fasta OUT.fasta SUPPORT.tsv MIN_TOOLS MIN_OVERLAP

IN.fasta holds the candidates of all plasmid predictors with headers
"tool|contig". A candidate is supported by a tool when that tool's candidates
cover >= MIN_OVERLAP of its length (nucmer alignments >= 95 % identity; the
own tool always counts). Candidates with support >= MIN_TOOLS go to OUT.fasta;
SUPPORT.tsv lists every candidate with its supporting tools.
Needs nucmer + show-coords (MUMmer) on PATH.
"""
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

MIN_IDY = 95.0


def read_fasta(path):
    recs, cur = [], None
    with open(path) as f:
        for ln in f:
            if ln.startswith(">"):
                cur = [ln[1:].strip().split()[0], []]
                recs.append(cur)
            elif cur:
                cur[1].append(ln.strip())
    return [(h, "".join(s)) for h, s in recs]


def covered(intervals):
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted(intervals):
        if cur_e is None or s > cur_e + 1:
            if cur_e is not None:
                total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s + 1
    return total


def main():
    inp, outp, supp, min_tools, min_ov = sys.argv[1:6]
    min_tools, min_ov = int(min_tools), float(min_ov)
    recs = read_fasta(inp)
    length = {h: len(s) for h, s in recs}
    tool = {h: h.split("|", 1)[0] for h, _ in recs}

    # intervals[query][other_tool] = [(start, end), ...] on the query
    intervals = defaultdict(lambda: defaultdict(list))
    with tempfile.TemporaryDirectory() as td:
        pfx = os.path.join(td, "all")
        # MUMmer 3 nucmer has no thread option (MUMmer 4: --threads)
        subprocess.run(["nucmer", "--maxmatch", "-p", pfx, inp, inp],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        out = subprocess.run(["show-coords", "-rclTH", pfx + ".delta"],
                             check=True, capture_output=True, text=True).stdout
    for ln in out.splitlines():
        p = ln.split("\t")
        if len(p) < 13:
            continue
        s1, e1, s2, e2 = (int(x) for x in p[0:4])
        idy, r, q = float(p[6]), p[11], p[12]
        if r == q or idy < MIN_IDY or tool.get(r) == tool.get(q):
            continue
        intervals[r][tool[q]].append((min(s1, e1), max(s1, e1)))
        intervals[q][tool[r]].append((min(s2, e2), max(s2, e2)))

    kept = 0
    with open(outp, "w") as fo, open(supp, "w") as fs:
        fs.write("candidate\ttool\tlength\tsupport\tsupporting_tools\tkept\n")
        for h, sq in recs:
            tools = {tool[h]}
            for t, iv in intervals[h].items():
                if covered(iv) >= min_ov * length[h]:
                    tools.add(t)
            ok = len(tools) >= min_tools
            kept += ok
            fs.write(f"{h}\t{tool[h]}\t{length[h]}\t{len(tools)}\t{','.join(sorted(tools))}\t{'yes' if ok else 'no'}\n")
            if ok:
                fo.write(f">{h}\n{sq}\n")
    print(f"[consensus] {kept}/{len(recs)} candidates supported by >= {min_tools} tool(s) "
          f"(overlap >= {min_ov:g} of length)")


if __name__ == "__main__":
    main()
