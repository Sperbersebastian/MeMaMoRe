#!/usr/bin/env python3
"""Fake bioinformatics tools for MeMaMoRe stub runs.

Every tool name in tests/stub/bin is a symlink to this file; we dispatch on
the name we were called as. Handlers write correctly named, plausible output
files so the real module scripts, .done bookkeeping and the GUI can run end to
end in seconds without envs or databases. They never read real data.
"""
import gzip
import json
import os
import random
import sys
from pathlib import Path

TOOL = os.path.basename(sys.argv[0])
ARGS = sys.argv[1:]

# Fixed toy community: (name, length, coverage, kind)
CONTIGS = [
    ("NODE_1_length_98000_cov_40.1", 98000, 40.1, "plasmid"),
    ("NODE_2_length_48500_cov_38.7", 48500, 38.7, "virus"),
    ("NODE_3_length_31000_cov_22.4", 31000, 22.4, "chromosome"),
    ("NODE_4_length_27000_cov_21.9", 27000, 21.9, "chromosome"),
    ("NODE_5_length_12000_cov_9.3", 12000, 9.3, "chromosome"),
    ("NODE_6_length_5400_cov_41.0", 5400, 41.0, "virus"),
    ("NODE_7_length_640_cov_3.0", 640, 3.0, "chromosome"),
]


def log(msg):
    sys.stderr.write(f"[stub:{TOOL}] {msg}\n")


def opt(*names, default=None):
    """Value following the first of `names` (also --name=value)."""
    for i, a in enumerate(ARGS):
        for n in names:
            if a == n and i + 1 < len(ARGS):
                return ARGS[i + 1]
            if n.startswith("--") and a.startswith(n + "="):
                return a.split("=", 1)[1]
    return default


def has(*names):
    return any(a in names for a in ARGS)


def positional(skip_values_for=()):
    out, skip = [], False
    for a in ARGS:
        if skip:
            skip = False
            continue
        if a in skip_values_for:
            skip = True
            continue
        if not a.startswith("-"):
            out.append(a)
    return out


def mkdir(p):
    Path(p).mkdir(parents=True, exist_ok=True)
    return Path(p)


def write(p, text):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)
    return p


def seq(name, length):
    rnd = random.Random(name)
    return "".join(rnd.choice("ACGT") for _ in range(length))


def write_fasta(p, records):
    lines = []
    for name, s in records:
        lines.append(f">{name}")
        lines.extend(s[i:i + 80] for i in range(0, len(s), 80))
    write(p, "\n".join(lines) + "\n")


def read_fasta_ids(p):
    ids = []
    try:
        with open(p) as f:
            for line in f:
                if line.startswith(">"):
                    ids.append(line[1:].split()[0])
    except OSError:
        pass
    return ids


def fasta_lengths(p):
    out, cur, n = {}, None, 0
    with open(p) as f:
        for line in f:
            if line.startswith(">"):
                if cur:
                    out[cur] = n
                cur, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if cur:
        out[cur] = n
    return out


def contig_info(cid):
    for name, length, cov, kind in CONTIGS:
        if cid == name or cid.startswith(name):
            return length, cov, kind
    return 1000, 5.0, "chromosome"


def copy_reads(src, dst):
    opener = gzip.open if str(src).endswith(".gz") else open
    with opener(src, "rb") as fi, gzip.open(dst, "wb") as fo:
        fo.write(fi.read())


def count_reads(p):
    opener = gzip.open if str(p).endswith(".gz") else open
    with opener(p, "rt") as f:
        return sum(1 for _ in f) // 4


# ---------------------------------------------------------------- QC
def fastp():
    i1, i2 = opt("-i", "--in1"), opt("-I", "--in2")
    o1, o2 = opt("-o", "--out1"), opt("-O", "--out2")
    copy_reads(i1, o1)
    n = count_reads(i1)
    if i2 and o2:
        copy_reads(i2, o2)
        n *= 2
    summ = {"total_reads": n, "total_bases": n * 150, "q30_rate": 0.93, "gc_content": 0.48}
    after = dict(summ, total_reads=int(n * 0.97), total_bases=int(n * 0.97) * 148)
    js = {
        "summary": {"before_filtering": summ, "after_filtering": after},
        "duplication": {"rate": 0.012},
        "adapter_cutting": {"adapter_trimmed_reads": n // 20},
        "filtering_result": {"passed_filter_reads": after["total_reads"],
                             "low_quality_reads": n - after["total_reads"], "too_short_reads": 0},
    }
    write(opt("-j", "--json"), json.dumps(js, indent=1))
    write(opt("-h", "--html"), "<html><body>stub fastp</body></html>\n")


def fastqc():
    out = mkdir(opt("-o", "--outdir"))
    for f in positional(skip_values_for=("-o", "--outdir", "-t", "--threads")):
        base = os.path.basename(f).replace(".fastq.gz", "").replace(".fq.gz", "")
        write(out / f"{base}_fastqc.html", "<html>stub fastqc</html>\n")
        write(out / f"{base}_fastqc.zip", "stub\n")


def multiqc():
    out = mkdir(opt("-o", "--outdir"))
    write(out / "multiqc_report.html", "<html>stub multiqc</html>\n")
    mkdir(out / "multiqc_data")


# ---------------------------------------------------------------- assembly
def spades():
    out = mkdir(opt("-o"))
    recs = [(n, seq(n, ln)) for n, ln, _, _ in CONTIGS]
    write_fasta(out / "contigs.fasta", recs)
    write_fasta(out / "scaffolds.fasta", recs)
    write(out / "spades.log", "stub spades\n")


def bwa_mem2():
    if ARGS and ARGS[0] == "index":
        ref = ARGS[-1]
        for ext in (".0123", ".amb", ".ann", ".pac", ".bwt.2bit.64"):
            write(ref + ext, "stub\n")
    elif ARGS and ARGS[0] == "mem":
        files = positional(skip_values_for=("-t", "-R", "-K"))[1:]
        sys.stdout.write(f"#STUBSAM\tref={os.path.abspath(files[0])}\n")


def bam_ref(bam):
    with open(bam) as f:
        first = f.readline().strip()
    return first.split("ref=", 1)[1] if "ref=" in first else None


def samtools():
    sub = ARGS[0] if ARGS else ""
    if sub == "sort":
        data = sys.stdin.read()
        write(opt("-o"), data or "#STUBSAM\n")
    elif sub == "index":
        write(positional()[-1] + ".bai", "stub\n")
    elif sub == "idxstats":
        ref = bam_ref(positional()[-1])
        for cid, ln in fasta_lengths(ref).items():
            _, cov, _ = contig_info(cid)
            print(f"{cid}\t{ln}\t{max(1, int(ln * cov / 150))}\t0")
        print("*\t0\t0\t120")
    elif sub == "view":
        if has("-c"):
            ref = bam_ref(positional()[-1])
            print(sum(int(ln * contig_info(c)[1] / 150) for c, ln in fasta_lengths(ref).items()) + 120)
        else:
            sys.stdout.write(sys.stdin.read())
    elif sub == "faidx":
        fa = positional()[1] if len(positional()) > 1 else positional()[0]
        with open(fa + ".fai", "w") as fo:
            for cid, ln in fasta_lengths(fa).items():
                fo.write(f"{cid}\t{ln}\t0\t80\t81\n")
    elif sub in ("flagstat", "stats"):
        print("1000 + 0 in total (QC-passed reads + QC-failed reads)")
    else:
        log(f"unhandled samtools {sub}; no-op")


def coverm():
    ref = opt("--reference", "-r")
    methods = []
    if "--methods" in ARGS or "-m" in ARGS:
        i = ARGS.index("--methods") if "--methods" in ARGS else ARGS.index("-m")
        for a in ARGS[i + 1:]:
            if a.startswith("-"):
                break
            methods.append(a)
    methods = methods or ["mean"]
    label = {"covered_fraction": "Covered Fraction", "mean": "Mean", "variance": "Variance",
             "length": "Length", "trimmed_mean": "Trimmed Mean", "count": "Read Count",
             "relative_abundance": "Relative Abundance (%)", "rpkm": "RPKM", "tpm": "TPM"}
    sample = "stub_R1.fastq.gz"
    rows = ["\t".join(["Contig"] + [f"{sample} {label.get(m, m)}" for m in methods])]
    ids = fasta_lengths(ref) if ref and os.path.exists(ref) else {}
    for cid, ln in ids.items():
        _, cov, _ = contig_info(cid)
        vals = {"covered_fraction": 0.99 if cov > 5 else 0.6, "mean": cov, "variance": cov / 4,
                "length": ln, "trimmed_mean": cov, "count": int(ln * cov / 150),
                "relative_abundance": round(100 / max(len(ids), 1), 3), "rpkm": cov * 10, "tpm": cov * 100}
        rows.append("\t".join([cid] + [str(vals.get(m, 0)) for m in methods]))
    out = opt("-o", "--output-file")
    text = "\n".join(rows) + "\n"
    if out:
        write(out, text)
    else:
        sys.stdout.write(text)


# ---------------------------------------------------------------- binning
# Toy binning result shared by all fake binners
BIN_PLAN = {"bin.1": ["NODE_3", "NODE_4"], "bin.2": ["NODE_5", "NODE_1"]}


def contigs_by_bin(fasta):
    """{bin: [(id, seq)]} for the contigs of `fasta` that BIN_PLAN assigns."""
    recs, cur, buf = [], None, []
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                if cur:
                    recs.append((cur, "".join(buf)))
                cur, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if cur:
        recs.append((cur, "".join(buf)))
    out = {}
    for b, prefixes in BIN_PLAN.items():
        for cid, sq in recs:
            if any(cid == p or cid.startswith(p + "_") for p in prefixes):
                out.setdefault(b, []).append((cid, sq))
    return out


def write_bins(fasta, outdir, prefix="", ext="fa"):
    mkdir(outdir)
    for b, recs in contigs_by_bin(fasta).items():
        write_fasta(Path(outdir) / f"{prefix}{b}.{ext}", recs)


def jgi_summarize_bam_contig_depths():
    bam = positional(skip_values_for=("--outputDepth", "--minMapQual", "--percentIdentity",
                                      "--referenceFasta"))[-1]
    ref = bam_ref(bam)
    name = os.path.basename(bam)
    rows = [f"contigName\tcontigLen\ttotalAvgDepth\t{name}\t{name}-var"]
    for cid, ln in fasta_lengths(ref).items():
        _, cov, _ = contig_info(cid)
        rows.append(f"{cid}\t{ln}\t{cov}\t{cov}\t{cov / 4}")
    write(opt("--outputDepth"), "\n".join(rows) + "\n")


def seqkit():
    if ARGS[:1] == ["seq"]:
        minlen = int(opt("-m", "--min-len", default="0"))
        fa = positional(skip_values_for=("-m", "--min-len", "-o"))[-1]
        recs, cur, buf = [], None, []
        with open(fa) as f:
            for line in f:
                if line.startswith(">"):
                    if cur:
                        recs.append((cur, "".join(buf)))
                    cur, buf = line[1:].rstrip("\n"), []
                else:
                    buf.append(line.strip())
        if cur:
            recs.append((cur, "".join(buf)))
        for name, sq in recs:
            if len(sq) >= minlen:
                sys.stdout.write(f">{name}\n{sq}\n")
    else:
        log(f"unhandled seqkit {ARGS[:1]}; no-op")


def metabat2():
    write_bins(opt("-i", "--inFile"), os.path.dirname(opt("-o", "--outFile")),
               prefix=os.path.basename(opt("-o", "--outFile")) + ".")
    # metabat2 names bins <prefix>.<n>.fa; our plan already uses "bin.N"
    out = Path(os.path.dirname(opt("-o", "--outFile")))
    for f in out.glob("bin.bin.*.fa"):
        f.rename(out / f.name.replace("bin.bin.", "bin.", 1))


def hmmpress():
    hmm = positional()[-1]
    for ext in ("h3f", "h3i", "h3m", "h3p"):
        write(f"{hmm}.{ext}", "stub\n")


def snakemake():
    """Binny: read --configfile, write <outputdir>/bins/*.fa."""
    cfg = opt("--configfile")
    conf = {}
    with open(cfg) as f:
        for line in f:
            line = line.strip()
            if ":" in line:
                k, v = line.split(":", 1)
                conf[k.strip()] = v.strip().strip('"')
    write_bins(conf["assembly"], Path(conf["outputdir"]) / "bins", prefix="binny_")


def run_comebin():
    out = Path(opt("-o"))
    write_bins(opt("-a"), out / "comebin_res" / "comebin_res_bins", prefix="comebin_")


def prodigal():
    i = opt("-i")
    prot = []
    for n, cid in enumerate(read_fasta_ids(i), 1):
        prot.append((f"{cid}_1 # 1 # 300 # 1 # ID=1_{n}", "M" + "A" * 99))
    if opt("-a"):
        write_fasta(opt("-a"), prot)
    if opt("-d"):
        write_fasta(opt("-d"), [(n, "ATG" + "GCT" * 99) for n, _ in prot])
    if opt("-o"):
        write(opt("-o"), "# stub prodigal\n")
    if opt("-f") == "gff" and not opt("-o"):
        sys.stdout.write("##gff-version 3\n")


def hmmsearch():
    faa = positional(skip_values_for=("-o", "--tblout", "--domtblout", "--cpu", "-E", "-T"))[-1]
    rows = ["# stub hmmsearch tblout"]
    for n, pid in enumerate(read_fasta_ids(faa)):
        acc = f"TIGR{n % 120:05d}" if "tigr" in (opt("--tblout") or "").lower() else f"PF{n % 120:05d}"
        rows.append(" ".join([pid, "-", acc, acc] + ["1e-50", "200.0", "0.0"] * 2 + ["1"] * 6 + ["-"]))
    for key in ("--tblout", "--domtblout"):
        if opt(key):
            write(opt(key), "\n".join(rows) + "\n")
    if opt("-o"):
        write(opt("-o"), "# stub\n")


def rscript():
    """MAGScoT.R -i <contig_to_bin> --hmm <map> -o <prefix>"""
    inp, prefix = opt("-i"), opt("-o")
    bins = {}
    with open(inp) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0].startswith("metabat2|"):
                bins.setdefault(parts[0], []).append(parts[1])
    if not bins:  # fall back to whatever binner is present
        with open(inp) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    bins.setdefault(parts[0], []).append(parts[1])
    rows, stats = ["binnew\tcontig"], ["bin\tcompleteness\tcontamination\tscore\trefined_id"]
    for i, (b, contigs) in enumerate(sorted(bins.items()), 1):
        rid = f"MAGScoT_cleanbin_{i:06d}"
        rows += [f"{rid}\t{c}" for c in contigs]
        stats.append(f"{b}\t0.92\t0.02\t0.88\t{rid}")
    write(prefix + ".refined.contig_to_bin.out", "\n".join(rows) + "\n")
    write(prefix + ".refined.out", "\n".join(stats) + "\n")
    write(prefix + ".scores.out", "\n".join(stats) + "\n")


def checkm2():
    if ARGS[:1] == ["database"]:
        write(Path(opt("--path")) / "CheckM2_database" / "uniref100.KO.1.dmnd", "stub\n")
        return
    out = mkdir(opt("--output-directory", "-o"))
    ext = opt("-x", "--extension", default="fa")
    rows = ["Name\tCompleteness\tContamination\tCompleteness_Model_Used\tTranslation_Table_Used\t"
            "Coding_Density\tContig_N50\tAverage_Gene_Length\tGenome_Size\tGC_Content\tTotal_Coding_Sequences\tAdditional_Notes"]
    for i, f in enumerate(sorted(Path(opt("--input", "-i")).glob(f"*.{ext}"))):
        size = sum(fasta_lengths(f).values())
        rows.append(f"{f.stem}\t{92.5 - 7 * i:.2f}\t{1.1 + i:.2f}\tNeural Network (Specific Model)\t11\t0.88\t{size}\t310\t{size}\t0.5\t{size // 1000}\tNone")
    write(out / "quality_report.tsv", "\n".join(rows) + "\n")


def barrnap():
    ids = [ln[1:].split()[0] for ln in sys.stdin if ln.startswith(">")]
    sys.stdout.write("##gff-version 3\n")
    if ids:
        sys.stdout.write(f"{ids[0]}\tbarrnap:0.9\trRNA\t100\t1600\t1e-200\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA\n")


def trnascan():
    out = opt("-o")
    fa = positional(skip_values_for=("-o",))[-1]
    ids = read_fasta_ids(fa)
    rows = ["Sequence\ttRNA\tBegin\tEnd\ttRNA\tAnti\tIntron Bounds\tInf", "Name\ttRNA #\tBounds\tBounds\tType\tCodon\tBegin\tEnd\tScore"]
    for i, cid in enumerate(ids[:1]):
        rows.append(f"{cid}\t1\t2000\t2073\tAla\tTGC\t0\t0\t70.1")
    write(out, "\n".join(rows) + "\n")


def gtdbtk():
    if ARGS[:1] == ["check_install"]:
        return
    if ARGS[:1] == ["classify_wf"]:
        out = mkdir(opt("--out_dir"))
        ext = opt("--extension", "-x", default="fa")
        rows = ["user_genome\tclassification\tfastani_reference\tclosest_placement_reference\tnote"]
        taxa = ["d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli",
                "d__Bacteria;p__Bacillota;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus;s__"]
        for i, f in enumerate(sorted(Path(opt("--genome_dir")).glob(f"*.{ext}"))):
            rows.append(f"{f.stem}\t{taxa[i % 2]}\tN/A\tN/A\tstub")
        write(out / "gtdbtk.bac120.summary.tsv", "\n".join(rows) + "\n")
        return
    log(f"unhandled gtdbtk {ARGS[:1]}; no-op")


# ---------------------------------------------------------------- plasmids / viruses
def records(fasta):
    out, cur, buf = [], None, []
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                if cur:
                    out.append((cur, "".join(buf)))
                cur, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if cur:
        out.append((cur, "".join(buf)))
    return out


def kind_of(cid):
    # ids may carry tool prefixes ("vv|NODE_1_...") from union steps
    return contig_info(cid.split("|")[-1])[2]


def metaplasmidspades():
    out = mkdir(opt("-o"))
    name, ln, _, _ = CONTIGS[0]
    write_fasta(out / "contigs.fasta", [(name + "_component_0", seq(name, ln))])


def metaviralspades():
    out = mkdir(opt("-o"))
    recs = [(n + "_component_0", seq(n, ln)) for n, ln, _, k in CONTIGS if k == "virus"]
    write_fasta(out / "contigs.fasta", recs)


def viralverify():
    fa, out = opt("-f"), mkdir(opt("-o"))
    base = Path(fa).name.rsplit(".", 1)[0]
    pred = {"plasmid": "Plasmid", "virus": "Virus", "chromosome": "Chromosome"}
    rows = ["Contig name,Prediction,Length,Score,Pfam hits"]
    by = {}
    for cid, sq in records(fa):
        p = pred[kind_of(cid)]
        rows.append(f"{cid},{p},{len(sq)},{12.5 if p != 'Chromosome' else -9.1},stub_hit")
        by.setdefault(p.lower(), []).append((cid, sq))
    write(out / f"{base}_result_table.csv", "\n".join(rows) + "\n")
    pdir = mkdir(out / "Prediction_results_fasta")
    for k in ("plasmid", "virus", "chromosome"):
        write_fasta(pdir / f"{base}_{k}.fasta", by.get(k, []))


def genomad():
    sub = ARGS[0] if ARGS else ""
    if sub == "download-database":
        dest = Path(positional()[-1]) / "genomad_db"
        write(dest / "version.txt", "1.9\n")
        return
    if sub != "end-to-end":
        log(f"unhandled genomad {sub}; no-op")
        return
    fa, out = positional(skip_values_for=("--threads", "-t", "--splits"))[1:3]
    base = Path(fa).name.rsplit(".", 1)[0]
    summ = mkdir(Path(out) / f"{base}_summary")
    prow = ["seq_name\tlength\ttopology\tcoordinates\tn_genes\tgenetic_code\tplasmid_score\tfdr\tn_hallmarks\tmarker_enrichment\tconjugation_genes\tamr_genes"]
    vrow = ["seq_name\tlength\ttopology\tcoordinates\tn_genes\tgenetic_code\tvirus_score\tfdr\tn_hallmarks\tmarker_enrichment\ttaxonomy"]
    pl, vi = [], []
    for cid, sq in records(fa):
        k = kind_of(cid)
        if k == "plasmid":
            prow.append(f"{cid}\t{len(sq)}\tNo terminal repeats\tNA\t{len(sq) // 1000}\t11\t0.9912\tNA\t3\t12.4\ttraA;traB\tNA")
            pl.append((cid, sq))
        elif k == "virus":
            vrow.append(f"{cid}\t{len(sq)}\tDTR\tNA\t{len(sq) // 1000}\t11\t0.9801\tNA\t5\t20.1\tViruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;")
            vi.append((cid, sq))
    write(summ / f"{base}_plasmid_summary.tsv", "\n".join(prow) + "\n")
    write(summ / f"{base}_virus_summary.tsv", "\n".join(vrow) + "\n")
    write_fasta(summ / f"{base}_plasmid.fna", pl)
    write_fasta(summ / f"{base}_virus.fna", vi)
    write(summ / f"{base}_summary.json", "{}\n")


MOB_COLS = ("sample_id\tnum_contigs\tsize\tgc\tmd5\trep_type(s)\trep_type_accession(s)\trelaxase_type(s)\t"
            "relaxase_type_accession(s)\tmpf_type\tmpf_type_accession(s)\torit_type(s)\torit_accession(s)\t"
            "predicted_mobility\tmash_nearest_neighbor\tmash_neighbor_distance\tmash_neighbor_identification\t"
            "primary_cluster_id\tsecondary_cluster_id\tpredicted_host_range_overall_rank\tpredicted_host_range_overall_name")


def mob_row(name, n, size):
    return (f"{name}\t{n}\t{size}\t50.1\tstubmd5\tIncFII\t000123__CP000001\tMOBF\tNC_000001\tMPF_F\tNC_000002\t"
            f"MOBF\tNC_000003\tconjugative\tCP000001\t0.01\tEscherichia coli\tAA001\tAA001\tfamily\tEnterobacteriaceae")


def mob_recon():
    fa, out = opt("-i", "--infile"), mkdir(opt("-o", "--outdir"))
    recs = records(fa)
    pl = [(c, s_) for c, s_ in recs if kind_of(c) == "plasmid"]
    ch = [(c, s_) for c, s_ in recs if kind_of(c) != "plasmid"]
    write_fasta(out / "plasmid_AA001.fasta", pl)
    write_fasta(out / "chromosome.fasta", ch)
    rep = ["sample_id\tmolecule_type\tprimary_cluster_id\tcontig_id\tsize"]
    rep += [f"stub\tplasmid\tAA001\t{c}\t{len(s_)}" for c, s_ in pl]
    rep += [f"stub\tchromosome\t-\t{c}\t{len(s_)}" for c, s_ in ch]
    write(out / "contig_report.txt", "\n".join(rep) + "\n")
    write(out / "mobtyper_results.txt", MOB_COLS + "\n" + mob_row("stub:AA001", len(pl), sum(len(x) for _, x in pl)) + "\n")


def mob_typer():
    fa, outfile = opt("-i", "--infile"), opt("-o", "--out_file")
    rows = [MOB_COLS] + [mob_row(c, 1, len(s_)) for c, s_ in records(fa)]
    write(outfile, "\n".join(rows) + "\n")


def virus_records(fasta):
    return [(c, s_) for c, s_ in records(fasta) if kind_of(c) == "virus"]


def virsorter():
    if ARGS[:1] == ["setup"]:
        write(Path(opt("-d", "--db-dir")) / "Done_all_setup", "stub\n")
        return
    fa, out = opt("-i", "--seqfile"), mkdir(opt("-w", "--working-dir"))
    vir = virus_records(fa)
    write_fasta(out / "final-viral-combined.fa", [(f"{c}||full", s_) for c, s_ in vir])
    rows = ["seqname\tdsDNAphage\tNCLDV\tRNA\tssDNA\tlavidaviridae\tmax_score\tmax_score_group\tlength\thallmark\tviral\tcellular"]
    rows += [f"{c}||full\t0.98\t0.1\t0.0\t0.2\t0.0\t0.98\tdsDNAphage\t{len(s_)}\t4\t60.0\t2.0" for c, s_ in vir]
    write(out / "final-viral-score.tsv", "\n".join(rows) + "\n")


def vibrant():
    fa, out = opt("-i"), opt("-folder")
    base = Path(fa).name.rsplit(".", 1)[0]
    d = mkdir(Path(out) / f"VIBRANT_{base}" / f"VIBRANT_phages_{base}")
    write_fasta(d / f"{base}.phages_combined.fna", virus_records(fa))


def checkv():
    if ARGS[:1] != ["end_to_end"]:
        log(f"unhandled checkv {ARGS[:1]}; no-op")
        return
    fa, out = positional(skip_values_for=("-t", "-d"))[1:3]
    out = mkdir(out)
    rows = ["contig_id\tcontig_length\tprovirus\tproviral_length\tgene_count\tviral_genes\thost_genes\t"
            "checkv_quality\tmiuvig_quality\tcompleteness\tcompleteness_method\tcontamination\tkmer_freq\twarnings"]
    for c, s_ in records(fa):
        q = "Complete" if len(s_) > 40000 else "High-quality"
        rows.append(f"{c}\t{len(s_)}\tNo\tNA\t{len(s_) // 1000}\t{len(s_) // 1500}\t0\t{q}\tHigh-quality\t100.0\tAAI-based (high-confidence)\t0.0\t1.0\t")
    write(out / "quality_summary.tsv", "\n".join(rows) + "\n")
    write(out / "completeness.tsv", "contig_id\tcompleteness\n")
    write(out / "contamination.tsv", "contig_id\tcontamination\n")


# ---------------------------------------------------------------- ARGs
def bwa():
    if ARGS and ARGS[0] == "index":
        ref = ARGS[-1]
        for ext in (".amb", ".ann", ".bwt", ".pac", ".sa"):
            write(ref + ext, "stub\n")
    elif ARGS and ARGS[0] == "mem":
        files = positional(skip_values_for=("-t", "-R", "-K"))[1:]
        sys.stdout.write(f"#STUBSAM\tref={os.path.abspath(files[0])}\n")


ARG_HITS = [("TEM-1", "beta-lactam", "ARO:3000873"), ("tetA", "tetracycline", "ARO:3000165")]


def deeparg():
    if ARGS[:1] == ["download_data"]:
        write(Path(opt("-o")) / "model" / "v2" / "metadata.pkl", "stub\n")
        return
    fa, prefix = opt("-i", "--input-file"), opt("-o", "--output-file")
    ids = read_fasta_ids(fa)
    hdr = ("#ARG\tquery-start\tquery-end\tread_id\tpredicted_ARG-class\tbest-hit\tprobability\t"
           "identity\talignment-length\talignment-bitscore\talignment-evalue\tcounts")
    rows = [hdr]
    for (gene, cls, _), cid in zip(ARG_HITS, ids):
        rows.append(f"{gene}\t10\t870\t{cid}\t{cls}\t{gene}|stub\t0.99\t98.5\t286\t560.0\t1e-150\t1")
    write(prefix + ".mapping.ARG", "\n".join(rows) + "\n")
    write(prefix + ".mapping.potential.ARG", hdr + "\n")


def rgi():
    sub = ARGS[0] if ARGS else ""
    if sub == "load":
        mkdir("localDB")
        write("localDB/card.json", "{}\n")
        return
    if sub != "main":
        log(f"unhandled rgi {sub}; no-op")
        return
    fa, prefix = opt("-i", "--input_sequence"), opt("-o", "--output_file")
    ids = read_fasta_ids(fa)
    hdr = ("ORF_ID\tContig\tStart\tStop\tOrientation\tCut_Off\tPass_Bitscore\tBest_Hit_Bitscore\t"
           "Best_Hit_ARO\tBest_Identities\tARO\tModel_type\tSNPs_in_Best_Hit_ARO\tOther_SNPs\t"
           "Drug Class\tResistance Mechanism\tAMR Gene Family\tPredicted_DNA\tPredicted_Protein\t"
           "CARD_Protein_Sequence\tPercentage Length of Reference Sequence\tID\tModel_ID\tNudged\tNote")
    rows = [hdr]
    for (gene, cls, aro), cid in zip(ARG_HITS, ids):
        rows.append(f"{cid}_1 # 1 # 861 # 1\t{cid}\t1\t861\t+\tStrict\t500\t560.0\t{gene}\t98.5\t"
                    f"{aro.split(':')[1]}\tprotein homolog model\tn/a\tn/a\t{cls} antibiotic\t"
                    f"antibiotic inactivation\tstub family\t\t\t\t100.0\tgnl|BL_ORD_ID|1\t1\t\t")
    write(prefix + ".txt", "\n".join(rows) + "\n")
    write(prefix + ".json", "{}\n")


HANDLERS = {
    "fastp": fastp, "fastqc": fastqc, "multiqc": multiqc,
    "metaspades.py": spades, "spades.py": spades,
    "bwa-mem2": bwa_mem2, "samtools": samtools, "coverm": coverm,
    "jgi_summarize_bam_contig_depths": jgi_summarize_bam_contig_depths, "seqkit": seqkit,
    "metabat2": metabat2, "hmmpress": hmmpress, "snakemake": snakemake,
    "run_comebin.sh": run_comebin, "prodigal": prodigal, "hmmsearch": hmmsearch,
    "Rscript": rscript, "checkm2": checkm2, "barrnap": barrnap, "tRNAscan-SE": trnascan,
    "gtdbtk": gtdbtk,
    "metaplasmidspades.py": metaplasmidspades, "metaviralspades.py": metaviralspades,
    "viralverify": viralverify, "genomad": genomad, "mob_recon": mob_recon, "mob_typer": mob_typer,
    "virsorter": virsorter, "VIBRANT_run.py": vibrant, "checkv": checkv,
    "bwa": bwa, "deeparg": deeparg, "rgi": rgi,
}

if __name__ == "__main__":
    if TOOL not in HANDLERS:
        log(f"no stub for '{TOOL}'")
        sys.exit(127)
    log(" ".join(ARGS)[:300])
    HANDLERS[TOOL]()
