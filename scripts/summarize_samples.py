#!/usr/bin/env python3
"""Multi-sample summary: collect per-sample results into comparison tables.

Usage: scripts/summarize_samples.py ROOT [--samples s1,s2,...] [--out DIR]
       (default samples: all in SRA/reads/manifest.tsv; default out: ROOT/SRA/summary)

Writes (tab-separated):
  samples_overview.tsv   one row per sample: reads, assembly, MAGs, plasmids, viruses, ARGs, virulence
  mags.tsv               every refined MAG: CheckM2, GUNC pass, quality tier, GTDB-Tk taxonomy
  viruses.tsv            every virus contig with CheckV quality (and iPHoP host genus if run)
  args_matrix.tsv        tool|gene x sample (hits over all contexts; RGI Loose and AMR++ <80% gene fraction excluded upstream)
  virulence_matrix.tsv   virulence gene x sample (VFDB via ABRicate)
  replicons_matrix.tsv   plasmid replicon x sample (PlasmidFinder via ABRicate)
Only the standard library is used; missing module outputs just leave cells empty.
"""
import argparse
import csv
import glob
import json
import os
from collections import defaultdict


def read_tsv(path):
    if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fasta_lengths(path):
    lens, n = [], None
    if not os.path.isfile(path):
        return lens
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if n is not None:
                    lens.append(n)
                n = 0
            elif n is not None:
                n += len(line.strip())
    if n is not None:
        lens.append(n)
    return lens


def n50(lens):
    total, acc = sum(lens), 0
    for x in sorted(lens, reverse=True):
        acc += x
        if acc >= total / 2:
            return x
    return 0


def count_fasta(path):
    return len(fasta_lengths(path)) if path and os.path.isfile(path) else ""


def mag_tier(comp, cont, gunc_pass=""):
    # MIMAG-style tiers from CheckM2 (completeness/contamination in %); a MAG
    # that GUNC flags as chimeric is "low" whatever its CheckM2 values
    if gunc_pass == "False":
        return "low"
    if comp >= 90 and cont < 5:
        return "high"
    if comp >= 50 and cont < 10:
        return "medium"
    return "low"


def samples_from_manifest(root):
    rows = read_tsv(os.path.join(root, "SRA", "reads", "manifest.tsv"))
    return [r["sample"] for r in rows if r.get("sample")]


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def write_matrix(path, label, counts, samples):
    genes = sorted(counts)
    write_tsv(path, [label] + samples,
              [[g] + [counts[g].get(s, 0) for s in samples] for g in genes])


def main():
    ap = argparse.ArgumentParser(description="MeMaMoRe multi-sample summary")
    ap.add_argument("root")
    ap.add_argument("--samples", help="comma-separated (default: all in SRA/reads/manifest.tsv)")
    ap.add_argument("--out", help="output dir (default: ROOT/SRA/summary)")
    a = ap.parse_args()

    root = os.path.abspath(a.root)
    sra = os.path.join(root, "SRA")
    samples = [s for s in (a.samples.split(",") if a.samples else samples_from_manifest(root)) if s]
    if not samples:
        raise SystemExit("[summary] no samples (empty manifest and no --samples)")
    out = a.out or os.path.join(sra, "summary")
    os.makedirs(out, exist_ok=True)

    overview, mags, viruses = [], [], []
    arg_m, vir_m, rep_m = defaultdict(dict), defaultdict(dict), defaultdict(dict)

    for s in samples:
        row = {"sample": s}

        # --- QC (fastp) ---
        fj = glob.glob(os.path.join(sra, "qc", "fastp", s, "*.fastp.json"))
        if fj:
            with open(fj[0]) as f:
                d = json.load(f).get("summary", {})
            row["reads_raw"] = d.get("before_filtering", {}).get("total_reads", "")
            row["reads_qc"] = d.get("after_filtering", {}).get("total_reads", "")
            q30 = d.get("after_filtering", {}).get("q30_rate")
            row["q30_qc"] = round(q30, 4) if isinstance(q30, (int, float)) else ""

        # --- assembly ---
        lens = fasta_lengths(os.path.join(sra, "assemblies", "spades", s, "contigs.fasta"))
        if lens:
            row["contigs"] = len(lens)
            row["contigs_ge1kb"] = sum(1 for x in lens if x >= 1000)
            row["assembly_bp"] = sum(lens)
            row["n50"] = n50(lens)
            row["longest"] = max(lens)

        # --- MAGs: CheckM2 + GTDB-Tk ---
        q = {r["Name"]: r for r in read_tsv(os.path.join(sra, "binning", "checkm2", s, "quality_report.tsv"))}
        tax = {}
        for f in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv"):
            for r in read_tsv(os.path.join(sra, "binning", "gtdbtk", s, f)):
                tax[r["user_genome"]] = r.get("classification", "")
        gunc = {}
        for gf in glob.glob(os.path.join(sra, "binning", "gunc", s, "GUNC.*.maxCSS_level.tsv")):
            for r in read_tsv(gf):
                gunc[r.get("genome", "")] = r.get("pass.GUNC", "")
        bins = sorted(glob.glob(os.path.join(sra, "binning", "magscot", s, "bins", "*.fa")))
        tiers = defaultdict(int)
        for b in bins:
            name = os.path.basename(b)[:-3]
            blens = fasta_lengths(b)
            r = q.get(name, {})
            try:
                comp, cont = float(r["Completeness"]), float(r["Contamination"])
                tier = mag_tier(comp, cont, gunc.get(name, ""))
            except (KeyError, ValueError):
                comp = cont = ""
                tier = "unknown"
            tiers[tier] += 1
            mags.append([s, name, len(blens), sum(blens), n50(blens), comp, cont, gunc.get(name, ""), tier,
                         tax.get(name, "")])
        if bins:
            row["mags"] = len(bins)
            row["mags_high"] = tiers["high"]
            row["mags_medium"] = tiers["medium"]

        # --- plasmids ---
        row["plasmids_derep"] = count_fasta(os.path.join(sra, "plasmids", s, "cluster", f"plasmids_derep.{s}.fasta"))

        # --- viruses: CheckV on the union set ---
        cv = read_tsv(os.path.join(sra, "viruses", s, "checkv", "quality_summary.tsv"))
        hosts = {}
        crispr_v = defaultdict(set)   # "<contig>" -> MAGs with a matching CRISPR spacer
        for r in read_tsv(os.path.join(sra, "summary", "crispr_hits.tsv")):
            if r.get("target_type") == "virus" and r.get("target", "").startswith(s + "|"):
                crispr_v[r["target"].split("|", 1)[1]].add(r.get("mag", ""))
        hp = os.path.join(sra, "viruses", s, "iphop", "host_genus.csv")
        if os.path.isfile(hp):
            with open(hp, newline="") as f:
                for r in csv.DictReader(f):
                    v, g = r.get("Virus", ""), r.get("Host genus", "")
                    if v and g and v not in hosts:  # first = best-scoring row
                        hosts[v] = f"{g} ({r.get('Confidence score', '')})"
        if cv:
            vq = defaultdict(int)
            for r in cv:
                qual = r.get("checkv_quality", "")
                vq[qual] += 1
                viruses.append([s, r.get("contig_id", ""), r.get("contig_length", ""), qual,
                                r.get("completeness", ""), r.get("provirus", ""),
                                hosts.get(r.get("contig_id", ""), ""),
                                ",".join(sorted(crispr_v.get(r.get("contig_id", ""), ())))])
            row["viruses"] = len(cv)
            row["viruses_complete_hq"] = vq["Complete"] + vq["High-quality"]

        # --- ARGs / virulence / replicons ---
        summ = os.path.join(sra, "args", s, "summary")
        args = read_tsv(os.path.join(summ, "args_summary.tsv"))
        if args:
            genes = set()
            for r in args:
                g = (r.get("gene") or "").strip()
                if g:
                    # tools name the same gene differently (CTX-M vs blaCTX-M-15_1): keep them apart
                    key = f"{r.get('source', '')}|{g}"
                    genes.add(key)
                    arg_m[key][s] = arg_m[key].get(s, 0) + 1
            row["arg_hits"] = len(args)
            row["arg_genes"] = len(genes)
        ctx = read_tsv(os.path.join(summ, "arg_context.tsv"))
        if ctx:
            row["arg_genes_ctx"] = len(ctx)
            row["arg_genes_plasmid"] = sum(1 for r in ctx if r.get("plasmid_cluster") or r.get("plasmid_tools"))
            row["arg_genes_mag"] = sum(1 for r in ctx if r.get("mag"))
            row["arg_genes_virus"] = sum(1 for r in ctx if r.get("virus_id"))
        vir = read_tsv(os.path.join(summ, "virulence_summary.tsv"))
        if vir:
            for r in vir:
                vir_m[r["gene"]][s] = vir_m[r["gene"]].get(s, 0) + 1
            row["virulence_genes"] = len({r["gene"] for r in vir})
        for r in read_tsv(os.path.join(summ, "plasmid_replicons.tsv")):
            rep_m[r["replicon"]][s] = rep_m[r["replicon"]].get(s, 0) + 1

        overview.append(row)

    cols = ["sample", "reads_raw", "reads_qc", "q30_qc", "contigs", "contigs_ge1kb", "assembly_bp",
            "n50", "longest", "mags", "mags_high", "mags_medium", "plasmids_derep", "viruses",
            "viruses_complete_hq", "arg_hits", "arg_genes", "arg_genes_ctx", "arg_genes_plasmid",
            "arg_genes_mag", "arg_genes_virus", "virulence_genes"]
    write_tsv(os.path.join(out, "samples_overview.tsv"), cols, [[r.get(c, "") for c in cols] for r in overview])
    write_tsv(os.path.join(out, "mags.tsv"),
              ["sample", "bin", "contigs", "size_bp", "n50", "completeness", "contamination", "gunc_pass", "quality",
               "gtdb_classification"],
              mags)
    write_tsv(os.path.join(out, "viruses.tsv"),
              ["sample", "contig", "length", "checkv_quality", "completeness", "provirus", "iphop_host_genus",
               "crispr_host_mags"],
              viruses)
    write_matrix(os.path.join(out, "args_matrix.tsv"), "tool|gene", arg_m, samples)
    write_matrix(os.path.join(out, "virulence_matrix.tsv"), "gene", vir_m, samples)
    write_matrix(os.path.join(out, "replicons_matrix.tsv"), "replicon", rep_m, samples)
    print(f"[summary] {len(samples)} sample(s), {len(mags)} MAGs, {len(viruses)} viruses, "
          f"{len(arg_m)} ARG genes, {len(vir_m)} virulence genes -> {out}")


if __name__ == "__main__":
    main()
