#!/usr/bin/env python3
"""Record tool, env, version and database version for a pipeline run.

Usage: collect_versions.py ROOT MODULE SAMPLE OUT.tsv [ENVS_JSON]

Versions come from the envs' conda-meta records (no tool is started), DB
versions from version files / directory names under refdata/. ENVS_JSON is the
output of `micromamba env list --json` (read from stdin when omitted / "-").
Writes OUT.tsv (module, sample, tool, env, package, version, build, database,
db_version) and appends the same rows, with a timestamp, to
ROOT/SRA/summary/versions_all.tsv; SRA/summary/versions.tsv keeps the newest
row per tool. Standard library only.
"""
import csv
import glob
import json
import os
import re
import sys
from datetime import datetime, timezone

# module -> [(tool, env, conda package)]
TOOLS = {
    "ingest": [("sra-tools", "env_ingest", "sra-tools"), ("pigz", "env_ingest", "pigz")],
    "qc": [("fastp", "env_qc_core", "fastp"), ("FastQC", "env_qc_core", "fastqc"),
           ("MultiQC", "env_qc_mqc", "multiqc"), ("Hostile", "hostile_env", "hostile")],
    "assembly": [("SPAdes", "env_assembly_core", "spades"), ("Pyrodigal", "env_assembly_core", "pyrodigal"),
                 ("bwa-mem2", "env_mapping_coverm", "bwa-mem2"), ("samtools", "env_mapping_coverm", "samtools"),
                 ("CoverM", "env_mapping_coverm", "coverm")],
    "binning": [("MetaBAT2", "env_binning", "metabat2"), ("Snakemake (Binny)", "env_binning", "snakemake"),
                ("COMEBin", "comebin_env", "comebin"), ("SemiBin2", "semibin_env", "semibin"),
                ("MAGScoT/R", "env_binning", "r-base"), ("HMMER", "env_binning", "hmmer"),
                ("GUNC", "gunc_env", "gunc"), ("CheckM2", "env_binning", "checkm2"),
                ("barrnap", "env_binning", "barrnap"), ("tRNAscan-SE", "env_binning", "trnascan-se"),
                ("GTDB-Tk", "env_binning", "gtdbtk")],
    "plasmids": [("metaPlasmidSPAdes", "env_assembly_core", "spades"), ("viralVerify", "viralverify_env", "viralverify"),
                 ("geNomad", "genomad_env", "genomad"), ("MOB-suite", "mobsuite_env", "mob_suite"),
                 ("Plasmer", "plasmer_env", "plasmer"), ("MUMmer", "stampede_env", "mummer"),
                 ("PLASMe/python", "plasme_env", "python"), ("HOTSPOT/python", "hotspot_env", "python")],
    "viruses": [("metaviralSPAdes", "env_assembly_core", "spades"), ("viralVerify", "viralverify_env", "viralverify"),
                ("geNomad", "genomad_env", "genomad"), ("VirSorter2", "virsorter2_env", "virsorter"),
                ("VIBRANT", "vibrant_env", "vibrant"), ("CheckV", "checkv_env", "checkv"),
                ("iPHoP", "iphop_env", "iphop")],
    "args": [("DeepARG", "deeparg_env", "deeparg"), ("RGI", "rgi_env", "rgi"), ("BWA (AMR++)", "amrplusplus_env", "bwa"),
             ("ABRicate", "abricate_env", "abricate"), ("AMRFinderPlus", "amrfinder_env", "ncbi-amrfinderplus"),
             ("argNorm", "amrfinder_env", "argnorm"), ("ISEScan", "mge_env", "isescan"),
             ("IntegronFinder", "mge_env", "integron_finder")],
}


def first_line(path):
    try:
        with open(path) as f:
            return f.readline().strip()
    except OSError:
        return ""


def newest(pattern):
    hits = sorted(glob.glob(pattern))
    return os.path.basename(hits[-1]) if hits else ""


def db_versions(root, envs):
    r = os.path.join(root, "refdata")
    card = ""
    try:
        with open(os.path.join(r, "card", "card.json")) as f:
            card = str(json.load(f).get("_version", ""))
    except (OSError, ValueError):
        pass
    abr = ""
    if "abricate_env" in envs:
        seqs = glob.glob(os.path.join(envs["abricate_env"], "db", "vfdb", "sequences"))
        if seqs:
            abr = "bundled DBs " + datetime.fromtimestamp(os.path.getmtime(seqs[0]), timezone.utc).strftime("%Y-%m-%d")
    amrf = ""
    if "amrfinder_env" in envs:
        amrf = newest(os.path.join(envs["amrfinder_env"], "share", "amrfinderplus", "data", "2*"))
    return {
        "geNomad": ("refdata/genomad/genomad_db", first_line(os.path.join(r, "genomad", "genomad_db", "version.txt"))),
        "GTDB-Tk": ("refdata/gtdbtk", newest(os.path.join(r, "gtdbtk", "release*"))),
        "CheckM2": ("refdata/checkm2", newest(os.path.join(r, "checkm2", "CheckM2_database", "*.dmnd"))),
        "GUNC": ("refdata/gunc", newest(os.path.join(r, "gunc", "gunc_db_*.dmnd"))),
        "CheckV": ("refdata/checkv", newest(os.path.join(r, "checkv", "checkv-db-*"))),
        "VirSorter2": ("refdata/virsorter2/db", "setup done" if os.path.exists(os.path.join(r, "virsorter2", "db", "Done_all_setup")) else ""),
        "VIBRANT": ("refdata/vibrant/db", "KEGG/Pfam32/VOG94" if os.path.exists(os.path.join(r, "vibrant", "db", "databases", "KEGG_profiles_prokaryotes.HMM.h3p")) else ""),
        "iPHoP": ("refdata/iphop", newest(os.path.join(r, "iphop", "*_rw"))),
        "DeepARG": ("refdata/deeparg", "model v2" if os.path.exists(os.path.join(r, "deeparg", "model", "v2", "model_LS.pkl")) else ""),
        "RGI": ("refdata/card", card and f"CARD {card}"),
        "BWA (AMR++)": ("refdata/amrplusplus", newest(os.path.join(r, "amrplusplus", "megares_v*"))),
        "AMRFinderPlus": ("env share/amrfinderplus/data", amrf),
        "ABRicate": ("env db/ (vfdb, resfinder, plasmidfinder)", abr),
        "PLASMe/python": ("refdata/plasme/DB", "present" if os.path.isdir(os.path.join(r, "plasme", "DB")) else ""),
        "Plasmer": ("refdata/plasmer", first_line(os.path.join(r, "plasmer", ".extracted"))),
        "MOB-suite": ("refdata/mobsuite_db", "present" if os.path.isdir(os.path.join(r, "mobsuite_db")) else ""),
        "Hostile": ("refdata/hostile", newest(os.path.join(r, "hostile", "*.1.bt2")).replace(".1.bt2", "")),
        "viralVerify": ("refdata/viralverify", "Pfam-A.hmm" if os.path.exists(os.path.join(r, "viralverify", "pfam", "Pfam-A.hmm")) else ""),
    }


def pkg_version(prefix, pkg):
    for meta in glob.glob(os.path.join(prefix, "conda-meta", f"{pkg}-*.json")):
        m = re.match(rf"^{re.escape(pkg)}-([^-]+)-([^-]+)\.json$", os.path.basename(meta))
        if m:
            return m.group(1), m.group(2)
    return "", ""


def main():
    root, module, sample, out = sys.argv[1:5]
    src = sys.argv[5] if len(sys.argv) > 5 else "-"
    data = json.load(sys.stdin if src == "-" else open(src))
    envs = {os.path.basename(p): p for p in data.get("envs", [])}
    dbs = db_versions(root, envs)
    modules = list(TOOLS) if module == "all" else [module]
    rows = []
    for mod in modules:
        for tool, env, pkg in TOOLS.get(mod, []):
            prefix = envs.get(env)
            ver, build = pkg_version(prefix, pkg) if prefix else ("", "")
            if not prefix:
                ver = "env missing"
            db, dbv = dbs.get(tool, ("", ""))
            rows.append([mod, sample, tool, env, pkg, ver, build, db, dbv])
    header = ["module", "sample", "tool", "env", "package", "version", "build", "database", "db_version"]
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)

    # project-wide history and newest row per (module, tool)
    sdir = os.path.join(root, "SRA", "summary")
    os.makedirs(sdir, exist_ok=True)
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    hist = os.path.join(sdir, "versions_all.tsv")
    new = not os.path.exists(hist)
    with open(hist, "a", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        if new:
            w.writerow(["timestamp"] + header)
        w.writerows([[ts] + r for r in rows])
    latest = {}
    with open(hist, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            latest[(r["module"], r["tool"])] = r
    with open(os.path.join(sdir, "versions.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["timestamp"] + header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(latest[k] for k in sorted(latest))
    print(f"[versions] {len(rows)} tools -> {out}")


if __name__ == "__main__":
    main()
