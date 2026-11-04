#!/usr/bin/env python3
"""Gene-centric ARG context table for one sample.

Usage: scripts/arg_context.py ROOT SAMPLE

Every ARG tool scans each gene (or, for ABRicate, each contig) once; this
script assigns the genomic context afterwards through the contig ID:

  MAG               SRA/binning/magscot/<s>/contigs_to_bin.tsv (refined bins)
  plasmid cluster   SRA/plasmids/<s>/cluster/*.clstr (Stampede, members "tool|contig")
  plasmid tools     SRA/plasmids/<s>/union/plasmids_raw.fasta headers ("tool|contig")
  virus             SRA/viruses/<s>/union/viruses_union.fasta; proviruses/fragments only
                    where the gene overlaps their coordinates (geNomad ID, VirSorter2
                    final-viral-boundary.tsv, VIBRANT prophage coordinates)

Writes SRA/args/<s>/summary/arg_context.tsv (one row per ARG gene) and
rewrites the 'context' column of args_summary.tsv (genes/contigs rows) to the
derived context (mag, plasmid, virus, unbinned; ';'-joined when several apply).

Hits per tool: DeepARG (genes.ffn), RGI Perfect/Strict (genes.faa; Loose only
with RGI_INCLUDE_LOOSE=1), AMRFinderPlus (AMR elements, genes.faa), ABRicate/
ResFinder (contigs; assigned to the gene with the largest overlap); 'aro' joins
the ARO terms of RGI and argNorm. AMR++ is read-based and has no gene context.
Standard library only.
"""
import csv
import os
import re
import sys
from collections import defaultdict


def read_tsv(path):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_tsv_commented(path):
    """TSV with leading '#' comment lines (argNorm output)."""
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return []
    with open(path, newline="") as f:
        lines = [ln for ln in f if not ln.startswith("# ")]
    return list(csv.DictReader(lines, delimiter="\t"))


def fasta_ids(path):
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        return [ln[1:].split()[0] for ln in f if ln.startswith(">")]


def load_genes(gff):
    """gene_id -> (contig, start, end); contig -> [(start, end, gene_id)]"""
    genes, by_contig = {}, defaultdict(list)
    if not os.path.isfile(gff):
        return genes, by_contig
    with open(gff) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "CDS":
                continue
            m = re.search(r"ID=([^;]+)", p[8])
            if not m:
                continue
            gid, start, end = m.group(1), int(p[3]), int(p[4])
            genes[gid] = (p[0], start, end)
            by_contig[p[0]].append((start, end, gid))
    return genes, by_contig


VIRUS_RE = re.compile(r"^(?P<contig>.+?)(?:\|provirus_(?P<s>\d+)_(?P<e>\d+)|\|\|.*|_fragment_\d+)?$")


def virus_coordinates(root, s):
    """virus_id -> (start, end) from the tools' own tables (VirSorter2 boundaries,
    VIBRANT prophage coordinates); geNomad proviruses carry them in the ID."""
    co = {}
    vdir = os.path.join(root, "SRA", "viruses", s)
    for r in read_tsv(os.path.join(vdir, "virsorter2", "final-viral-boundary.tsv")):
        try:
            co[r["seqname_new"]] = (int(r["trim_bp_start"]), int(r["trim_bp_end"]))
        except (KeyError, ValueError):
            pass
    for dirpath, _, files in os.walk(os.path.join(vdir, "vibrant")):
        for fn in files:
            if fn.startswith("VIBRANT_integrated_prophage_coordinates_") and fn.endswith(".tsv"):
                for r in read_tsv(os.path.join(dirpath, fn)):
                    try:
                        co[r["fragment"]] = (int(r["nucleotide start"]), int(r["nucleotide stop"]))
                    except (KeyError, ValueError):
                        pass
    return co


def load_viruses(fasta, coords):
    """contig -> [(start or None, end or None, virus_id)]; None = whole contig"""
    out = defaultdict(list)
    for vid in fasta_ids(fasta):
        m = VIRUS_RE.match(vid)
        c = m.group("contig")
        if m.group("s"):
            s, e = int(m.group("s")), int(m.group("e"))
        else:
            s, e = coords.get(vid, (None, None))
        out[c].append((s, e, vid))
    return out


def load_plasmids(root, s):
    cluster, tools = {}, defaultdict(set)
    for clstr in sorted(os.listdir(os.path.join(root, "SRA", "plasmids", s, "cluster"))
                        if os.path.isdir(os.path.join(root, "SRA", "plasmids", s, "cluster")) else []):
        if not clstr.endswith(".clstr"):
            continue
        cur = None
        with open(os.path.join(root, "SRA", "plasmids", s, "cluster", clstr)) as f:
            for ln in f:
                p = ln.rstrip("\n").split("\t")
                if p[0].startswith(">"):
                    cur = p[0][1:]
                    member = p[1] if len(p) > 1 else ""
                else:
                    member = p[0]
                if cur and "|" in member:
                    cluster.setdefault(member.split("|", 1)[1], cur)
    for h in fasta_ids(os.path.join(root, "SRA", "plasmids", s, "union", "plasmids_raw.fasta")):
        if "|" in h:
            tool, contig = h.split("|", 1)
            tools[contig].add(tool)
    return cluster, tools


def main():
    root, s = os.path.abspath(sys.argv[1]), sys.argv[2]
    args_dir = os.path.join(root, "SRA", "args", s)
    out_dir = os.path.join(args_dir, "summary")
    os.makedirs(out_dir, exist_ok=True)

    genes, by_contig = load_genes(os.path.join(root, "SRA", "assemblies", "genes", s, "genes.gff"))
    mag = {r["contig"]: r["bin"] for r in read_tsv(os.path.join(root, "SRA", "binning", "magscot", s, "contigs_to_bin.tsv"))}
    pcluster, ptools = load_plasmids(root, s)
    viruses = load_viruses(os.path.join(root, "SRA", "viruses", s, "union", "viruses_union.fasta"),
                           virus_coordinates(root, s))

    def contig_of(gid):
        return genes[gid][0] if gid in genes else re.sub(r"_\d+$", "", gid)

    def virus_of(gid, contig):
        hits = []
        for vs, ve, vid in viruses.get(contig, []):
            if vs is None or gid not in genes:
                hits.append(vid)
            else:
                _, gs, ge = genes[gid]
                if gs <= ve and ge >= vs:
                    hits.append(vid)
        return ";".join(sorted(set(hits)))

    # tool -> gene -> label ("name (class)")
    hits = defaultdict(dict)
    for r in read_tsv(os.path.join(args_dir, "deeparg", "genes", "genes.mapping.ARG")):
        gid = (r.get("read_id") or "").split()[0]
        if gid:
            hits["deeparg"][gid] = f"{r.get('#ARG', '')} ({r.get('predicted_ARG-class', '')})"
    loose = os.environ.get("RGI_INCLUDE_LOOSE", "0") == "1"
    for r in read_tsv(os.path.join(args_dir, "rgi", "genes", "genes.txt")):
        if r.get("Cut_Off") == "Loose" and not loose:
            continue
        gid = (r.get("ORF_ID") or "").split()[0]
        if gid:
            hits["rgi"][gid] = f"{r.get('Best_Hit_ARO', '')} ({r.get('Drug Class', '')})"
    for r in read_tsv(os.path.join(args_dir, "amrfinderplus", "genes.tsv")):
        if r.get("Type") == "AMR" and r.get("Protein id"):
            hits["amrfinderplus"][r["Protein id"]] = f"{r.get('Element symbol', '')} ({r.get('Class', '')})"

    # ARO terms: RGI reports them, argNorm adds them for DeepARG, AMRFinderPlus, ResFinder
    aro = defaultdict(set)

    def add_aro(gid, value):
        for a in (value or "").replace(";", ",").split(","):
            a = a.strip()
            if a.startswith("ARO:"):
                aro[gid].add(a)
            elif a.isdigit():
                aro[gid].add("ARO:" + a)
    for r in read_tsv(os.path.join(args_dir, "rgi", "genes", "genes.txt")):
        if r.get("Cut_Off") != "Loose" or loose:
            add_aro((r.get("ORF_ID") or "").split()[0], r.get("ARO"))
    for r in read_tsv_commented(os.path.join(args_dir, "deeparg", "genes", "genes.argnorm.tsv")):
        add_aro((r.get("read_id") or "").split()[0], r.get("ARO"))
    for r in read_tsv_commented(os.path.join(args_dir, "amrfinderplus", "genes.argnorm.tsv")):
        if r.get("Type") == "AMR":
            add_aro(r.get("Protein id", ""), r.get("ARO"))
    abr_aro = {(r.get("SEQUENCE"), r.get("START"), r.get("END")): r.get("ARO")
               for r in read_tsv_commented(os.path.join(args_dir, "abricate", "resfinder", "contigs", "contigs.argnorm.tsv"))}

    for r in read_tsv(os.path.join(args_dir, "abricate", "resfinder", "contigs", "contigs.tsv")):
        c, hs, he = r.get("SEQUENCE", ""), int(r.get("START") or 0), int(r.get("END") or 0)
        best, best_ov = None, 0
        for gs, ge, gid in by_contig.get(c, []):
            ov = min(ge, he) - max(gs, hs) + 1
            if ov > best_ov:
                best, best_ov = gid, ov
        if best:
            hits["abricate_resfinder"][best] = f"{r.get('GENE', '')} ({r.get('RESISTANCE', '')})"
            add_aro(best, abr_aro.get((r.get("SEQUENCE"), r.get("START"), r.get("END"))))

    tools = ["deeparg", "rgi", "amrfinderplus", "abricate_resfinder"]
    all_genes = sorted({g for t in tools for g in hits[t]},
                       key=lambda g: (contig_of(g), genes.get(g, ("", 0, 0))[1]))
    cols = ["gene_id", "contig", "start", "end", "mag", "plasmid_cluster", "plasmid_tools", "virus_id",
            "n_tools", "aro"] + tools + ["plasmid_contig_in_mag"]
    rows = []
    for g in all_genes:
        c = contig_of(g)
        _, gs, ge = genes.get(g, (c, "", ""))
        pc = pcluster.get(c, "")
        pt = ",".join(sorted(ptools.get(c, ())))
        m = mag.get(c, "")
        rows.append({
            "gene_id": g, "contig": c, "start": gs, "end": ge, "mag": m,
            "plasmid_cluster": pc, "plasmid_tools": pt, "virus_id": virus_of(g, c),
            "n_tools": sum(1 for t in tools if g in hits[t]),
            "aro": ",".join(sorted(aro.get(g, ()))),
            **{t: hits[t].get(g, "") for t in tools},
            "plasmid_contig_in_mag": "yes" if (m and (pc or pt)) else "",
        })
    with open(os.path.join(out_dir, "arg_context.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    # derived context for args_summary.tsv rows that point at a gene or contig
    summ = os.path.join(out_dir, "args_summary.tsv")
    srows = read_tsv(summ)
    if srows:
        def label(contig, gid=None):
            parts = []
            if mag.get(contig):
                parts.append("mag")
            if pcluster.get(contig) or ptools.get(contig):
                parts.append("plasmid")
            if virus_of(gid or contig, contig):
                parts.append("virus")
            return ";".join(parts) or "unbinned"
        for r in srows:
            if r.get("context") in ("genes", "contigs"):
                ref = (r.get("contig_or_gene") or "").split()[0]
                contig = contig_of(ref) if ref in genes else ref
                r["context"] = label(contig, ref if ref in genes else None)
        with open(summ, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(srows[0].keys()), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(srows)

    # ABRicate VFDB / PlasmidFinder tables: context from the contig as well
    for fname in ("virulence_summary.tsv", "plasmid_replicons.tsv"):
        p = os.path.join(out_dir, fname)
        trows = read_tsv(p)
        if not trows:
            continue
        for r in trows:
            if r.get("context") == "contigs":
                c = r.get("contig", "")
                parts = [x for x, ok in (("mag", mag.get(c)), ("plasmid", pcluster.get(c) or ptools.get(c)),
                                         ("virus", viruses.get(c))) if ok]
                r["context"] = ";".join(parts) or "unbinned"
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(trows[0].keys()), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(trows)

    n_pl = sum(1 for r in rows if r["plasmid_cluster"] or r["plasmid_tools"])
    print(f"[arg-context] {len(rows)} ARG genes ({n_pl} on plasmid contigs, "
          f"{sum(1 for r in rows if r['mag'])} in MAGs, {sum(1 for r in rows if r['virus_id'])} on viruses, "
          f"{sum(1 for r in rows if r['plasmid_contig_in_mag'])} plasmid-contig-in-MAG) -> {out_dir}/arg_context.tsv")


if __name__ == "__main__":
    main()
