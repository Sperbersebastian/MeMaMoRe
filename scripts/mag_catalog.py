#!/usr/bin/env python3
"""Build SRA/summary/mag_catalog.tsv from a dRep v4 run and CoverM abundances.

Usage: mag_catalog.py ROOT CATALOG_DIR SAMPLE [SAMPLE ...]

One row per representative genome (dRep winner): secondary cluster, members,
CheckM2 completeness/contamination, GUNC pass and GTDB-Tk taxonomy of the
representative, then per sample relative abundance (%) and mean coverage
from CoverM genome. Standard library only.
"""
import csv
import glob
import os
import sys
from collections import defaultdict


def read(path, delim=","):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter=delim))


def main():
    root, cat, samples = os.path.abspath(sys.argv[1]), sys.argv[2], sys.argv[3:]
    sra = os.path.join(root, "SRA")
    drep = os.path.join(cat, "drep", "data_tables")
    cdb = read(os.path.join(drep, "Cdb.csv"))
    winners = {r["genome"]: r for r in read(os.path.join(drep, "Wdb.csv"))}
    info = {r["genome"]: r for r in read(os.path.join(cat, "genomeInfo.csv"))}
    members = defaultdict(list)
    for r in cdb:
        members[r["secondary_cluster"]].append(r["genome"])

    def split(genome):                       # "<sample>__<bin>.fa"
        s, b = genome.split("__", 1)
        return s, b[:-3] if b.endswith(".fa") else b

    tax, gunc = {}, {}
    for s in set(samples):
        for f in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv"):
            for r in read(os.path.join(sra, "binning", "gtdbtk", s, f), "\t"):
                tax[(s, r["user_genome"])] = r.get("classification", "")
        for gf in glob.glob(os.path.join(sra, "binning", "gunc", s, "GUNC.*.maxCSS_level.tsv")):
            for r in read(gf, "\t"):
                gunc[(s, r.get("genome", ""))] = r.get("pass.GUNC", "")

    # CoverM genome: Genome, "<reads> Relative Abundance (%)", "<reads> Mean", "<reads> Covered Fraction"
    ab = defaultdict(dict)
    for s in samples:
        for r in read(os.path.join(cat, "abundance", f"{s}.tsv"), "\t"):
            g = r.get("Genome", "")
            for k, v in r.items():
                if k.endswith("Relative Abundance (%)"):
                    ab[g][f"{s}_rel_abundance_pct"] = v
                elif k.endswith(" Mean"):
                    ab[g][f"{s}_mean_coverage"] = v

    cols = ["catalog_genome", "secondary_cluster", "n_members", "members", "sample", "bin",
            "completeness", "contamination", "gunc_pass", "gtdb_classification"]
    for s in samples:
        cols += [f"{s}_rel_abundance_pct", f"{s}_mean_coverage"]
    rows = []
    for g, w in sorted(winners.items(), key=lambda kv: kv[1].get("cluster", "")):
        s, b = split(g)
        cl = w.get("cluster", "")
        mem = sorted(members.get(cl, [g]))
        row = {"catalog_genome": g[:-3] if g.endswith(".fa") else g, "secondary_cluster": cl,
               "n_members": len(mem), "members": ",".join(m[:-3] for m in mem), "sample": s, "bin": b,
               "completeness": info.get(g, {}).get("completeness", ""),
               "contamination": info.get(g, {}).get("contamination", ""),
               "gunc_pass": gunc.get((s, b), ""), "gtdb_classification": tax.get((s, b), "")}
        stem = g[:-3] if g.endswith(".fa") else g
        for k in cols[10:]:
            row[k] = ab.get(stem, {}).get(k, "")
        rows.append(row)
    out = os.path.join(sra, "summary", "mag_catalog.tsv")
    with open(out, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        wr.writeheader()
        wr.writerows(rows)
    print(f"[catalog] {len(rows)} catalog genome(s), {sum(len(v) for v in members.values())} MAGs -> {out}")


if __name__ == "__main__":
    main()
