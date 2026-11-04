#!/usr/bin/env python3
"""Plasmid host evidence table: HOTSPOT next to abundance-profile correlation.

Usage: plasmid_hosts.py ROOT CATALOG_DIR MIN_SAMPLES MIN_RHO MAX_P SAMPLE [SAMPLE ...]

For every global plasmid cluster representative (SRA/plasmids/_global/cluster/
plasmids_derep.global.fasta, headers "<sample>|<tool>|<contig>"):
  hotspot_*      HOTSPOT host lineage predicted in the plasmid's own sample
  corr_*         MAG catalog genome whose relative abundance across samples
                 correlates best (Spearman) with the plasmid's abundance (CoverM
                 TPM per sample, CATALOG_DIR/plasmid_abundance/<s>.tsv); reported
                 when rho >= MIN_RHO and permutation p <= MAX_P, and only with
                 >= MIN_SAMPLES samples (otherwise a note)
  crispr_*       filled by the CRISPR spacer step when present (crispr_hits.tsv)
Writes SRA/summary/plasmid_hosts.tsv. Standard library only.
"""
import csv
import os
import random
import sys


def read(path, delim="\t"):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter=delim))


def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sx = sum((a - mx) ** 2 for a in x) ** 0.5
    sy = sum((b - my) ** 2 for b in y) ** 0.5
    return sxy / (sx * sy) if sx and sy else 0.0


def spearman(x, y, n_perm=2000, seed=1):
    rx, ry = ranks(x), ranks(y)
    rho = pearson(rx, ry)
    rnd = random.Random(seed)
    hits, ry2 = 0, list(ry)
    for _ in range(n_perm):
        rnd.shuffle(ry2)
        if pearson(rx, ry2) >= rho - 1e-12:
            hits += 1
    return rho, (hits + 1) / (n_perm + 1)


def main():
    root, cat = os.path.abspath(sys.argv[1]), sys.argv[2]
    min_samples, min_rho, max_p = int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
    samples = sys.argv[6:]
    sra = os.path.join(root, "SRA")
    reps = []
    gfa = os.path.join(sra, "plasmids", "_global", "cluster", "plasmids_derep.global.fasta")
    if os.path.isfile(gfa):
        cur = None
        with open(gfa) as f:
            for ln in f:
                if ln.startswith(">"):
                    cur = [ln[1:].split()[0], 0]
                    reps.append(cur)
                elif cur:
                    cur[1] += len(ln.strip())

    # HOTSPOT per sample: Contig "<tool>|<contig>"
    hot = {}
    for s in samples:
        for r in read(os.path.join(sra, "plasmids", s, "hotspot", "host_lineage.tsv")):
            hot[(s, r.get("Contig", ""))] = r

    # abundance profiles
    pab = {}   # rep -> {sample: tpm}
    for s in samples:
        for r in read(os.path.join(cat, "plasmid_abundance", f"{s}.tsv")):
            c = r.get("Contig", "")
            for k, v in r.items():
                if k.endswith(" TPM"):
                    pab.setdefault(c, {})[s] = float(v or 0)
    mags = read(os.path.join(sra, "summary", "mag_catalog.tsv"))
    mab = {}
    for m in mags:
        mab[m["catalog_genome"]] = {s: float(m.get(f"{s}_rel_abundance_pct") or 0) for s in samples}
    tax = {m["catalog_genome"]: m.get("gtdb_classification", "") for m in mags}

    crispr = {}
    for r in read(os.path.join(sra, "summary", "crispr_hits.tsv")):
        if r.get("target_type") == "plasmid":
            crispr.setdefault(r["target"], set()).add(r["mag"])

    enough = len(samples) >= min_samples
    cols = ["plasmid", "sample", "length", "n_samples_detected", "hotspot_order", "hotspot_family", "hotspot_genus",
            "corr_host_genome", "corr_rho", "corr_p", "corr_host_taxonomy", "corr_note", "crispr_host_mags"]
    rows = []
    for rep, length in reps:
        s, rest = rep.split("|", 1) if "|" in rep else ("", rep)
        h = hot.get((s, rest), {})
        prof = pab.get(rep, {})
        x = [prof.get(t, 0.0) for t in samples]
        row = {"plasmid": rep, "sample": s, "length": length,
               "n_samples_detected": sum(1 for v in x if v > 0),
               "hotspot_order": h.get("Order", ""), "hotspot_family": h.get("Family", ""),
               "hotspot_genus": h.get("Genus", ""), "corr_host_genome": "", "corr_rho": "", "corr_p": "",
               "corr_host_taxonomy": "", "corr_note": "", "crispr_host_mags": ",".join(sorted(crispr.get(rep, ())))}
        if not enough:
            row["corr_note"] = f"needs >= {min_samples} samples (have {len(samples)})"
        elif not prof:
            row["corr_note"] = "no abundance profile"
        else:
            best = None
            for g, p in mab.items():
                y = [p.get(t, 0.0) for t in samples]
                if len(set(x)) < 2 or len(set(y)) < 2:
                    continue
                rho, pv = spearman(x, y)
                if best is None or rho > best[1]:
                    best = (g, rho, pv)
            if best and best[1] >= min_rho and best[2] <= max_p:
                row.update(corr_host_genome=best[0], corr_rho=f"{best[1]:.3f}", corr_p=f"{best[2]:.4f}",
                           corr_host_taxonomy=tax.get(best[0], ""))
            else:
                row["corr_note"] = "no MAG profile passes rho/p thresholds" if best else "constant profile"
        rows.append(row)
    out = os.path.join(sra, "summary", "plasmid_hosts.tsv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    n_corr = sum(1 for r in rows if r["corr_host_genome"])
    print(f"[hosts] {len(rows)} plasmid clusters, {n_corr} with a correlated MAG host"
          f"{'' if enough else f' (correlation needs >= {min_samples} samples)'} -> {out}")


if __name__ == "__main__":
    main()
