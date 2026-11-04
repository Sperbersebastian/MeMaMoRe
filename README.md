# **MeMaMoRe: Metagenomic MAG/Mobilome/Resistome Reconstruction**

**Version:** 0.9-dev
**Author:** Sebastian Sperber
**Institution:** University of Potsdam & Leibniz Institute for Agricultural Engineering and Bioeconomy (ATB)
**License:** MIT
**Repository:** [github.com/Sperbersebastian/MeMaMoRe](https://github.com/Sperbersebastian/MeMaMoRe)
**Contact:** [se.sperber@gmail.com](mailto:se.sperber@gmail.com)

---

## **Overview**

**MeMaMoRe** is a modular and reproducible framework for the reconstruction and characterization of **metagenome-assembled genomes (MAGs)**, **plasmids**, and **viruses** from shotgun metagenomic data.
It provides a complete, tool-integrated analysis chain with high flexibility and modular environment control.
The pipeline supports **local execution** and **HPC clusters**, using **Micromamba** for isolated environment management.

Core components include:

* **Quality Control and Trimming**
* **Assembly**
* **Binning and Refinement**
* **MAG Quality and Taxonomy**
* **Plasmid and Virus Detection**
* **ARG and Virulence Profiling**
* **(optional)** Simulation for test and benchmarking datasets

---

## **Key Features**

| Module           | Description                                          | Tools                                                                          |
| ---------------- | ---------------------------------------------------- | ------------------------------------------------------------------------------ |
| **QC**           | Read quality trimming and report generation          | Fastp, FastQC, MultiQC                                                         |
| **Assembly**     | Metagenomic and plasmid-targeted assembly            | metaSPAdes, metaPlasmidSPAdes                                                  |
| **Binning**      | Genome binning and refinement                        | MetaBAT2, Binny, COMEBin, MAGScoT, CheckM2                                     |
| **Taxonomy**     | Genome classification via GTDB                       | GTDB-Tk (r226)                                                                 |
| **Plasmids**     | Detection, clustering, typing and host prediction    | metaPlasmidSPAdes, ViralVerify, geNomad, PLASMe, MOB-suite, Stampede, HOTSPOT  |
| **Viruses**      | Viral prediction, dereplication, QC, host prediction | metaviralSPAdes, ViralVerify, geNomad, VirSorter2, VIBRANT, CheckV, iPHoP (optional) |
| **ARGs**         | Resistance genes, virulence factors, replicons       | DeepARG, RGI/CARD, AMR++/MEGARes, ABRicate (VFDB, ResFinder, PlasmidFinder)    |
| **Summary**      | Multi-sample comparison tables                       | `bin/main.sh summary` → `SRA/summary/`                                         |
| **Simulation** ⚠ | Synthetic reads and controlled variation *(partial)* | ART, custom `sim_env`                                                          |
| **Web GUI**      | Dashboard for sample management and visualization    | Flask, HTML/JS/CSS                                                             |

---

## **Directory Structure**

```
MeMaMoRe/
├── bin/
│   ├── main.sh              # Main CLI entry
│   ├── lib/                 # Shared helper functions (env mgmt, logging)
│   └── modules/
│       ├── ingest.sh        # Sample import / metadata management
│       ├── qc.sh            # Quality control orchestrator
│       ├── assembly.sh      # Assembly orchestrator
│       ├── binning.sh       # Binning/refinement orchestrator
│       ├── plasmids.sh      # Plasmid workflow orchestrator
│       ├── viruses.sh       # Virus workflow orchestrator
│       ├── args.sh          # ARG / virulence workflow orchestrator
│       └── sim.sh           # Simulation orchestrator
├── envs/                    # Micromamba environment creation scripts
│   ├── create_env_plasmids.sh
│   ├── create_env_binning.sh
│   ├── create_env_qc.sh
│   └── ...
├── external/                # External cloned tools (e.g. Binny, COMEBin, HOTSPOT)
├── gui/                     # Web Dashboard (Flask application)
├── memamore-gui.sh          # Launcher script for the Web GUI
├── modules/                 # Tool-specific executors (assembly, plasmids, binning)
├── refdata/                 # Databases (GTDB, PLASMe, geNomad, HOTSPOT, etc.)
├── SRA/                     # Main workspace for input/output per sample
│   ├── reads/
│   ├── assemblies/
│   ├── binning/
│   ├── plasmids/
│   ├── viruses/
│   ├── args/
│   ├── qc/
│   └── summary/             # multi-sample tables (bin/main.sh summary)
├── scripts/                 # helpers (params merge, test data, multi-sample summary)
├── tests/                   # GUI smoke test, stub pipeline (fake tools) for CI
└── config/                  # YAML configurations (params, testsets)
```

⚠ The **simulation** module is partially implemented.

---

## **Installation**

### 1. Prerequisites

* Linux (tested on Ubuntu 22.04, 25.x and CentOS 7); Windows via WSL2
* Micromamba (`curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C ~/micromamba bin/micromamba`)
* Git ≥ 2.30
* Bash ≥ 5.0

### 2. Clone the Repository

```bash
git clone https://github.com/Sperbersebastian/MeMaMoRe.git
cd MeMaMoRe
```

### 3. Create Environments

```bash
bin/main.sh env create all          # or one module: ingest|qc|assembly|binning|plasmids|viruses|args|...
bin/main.sh env list                # module -> env -> creator script
```

Missing environments are also created on first use. Reference databases live in
`refdata/` (override with `REFDATA_BASE`, `GTDBTK_DATA_PATH`, `VIBRANT_DB`,
`IPHOP_DB`); most are downloaded on first run. Rough needs for everything:
~300 GB disk (plus ~340 GB for the optional iPHoP DB), and ≥64 GB RAM for
GTDB-Tk (skipped with a message on smaller machines, `GTDBTK_FORCE=1` overrides).

---

## **Usage**

### General Syntax

```bash
bin/main.sh run <module> --sample <SAMPLE> [--force] [--force-env]
```

### Examples

**Run quality control:**

```bash
bin/main.sh run qc --sample test_sample
```

**Run assembly:**

```bash
bin/main.sh run assembly --sample test_sample
```

**Run plasmid analysis:**

```bash
bin/main.sh run plasmids --sample test_sample
```

**Global plasmid dereplication:**

```bash
bin/main.sh run plasmids --global-cluster
```

**Run binning workflow:**

```bash
bin/main.sh run binning --sample test_sample
```

**Viruses, ARGs/virulence, everything:**

```bash
bin/main.sh run viruses --sample test_sample
bin/main.sh run args --sample test_sample
bin/main.sh run all --sample test_sample
```

**Compare all samples** (tables in `SRA/summary/`, also shown on the GUI dashboard):

```bash
bin/main.sh summary                 # or: --samples s1,s2
```

---

## **Web GUI (Dashboard)**

MeMaMoRe includes a fully-featured Web GUI for easier management of samples, pipeline execution, and result visualization.

### Features
* **Dashboard:** Overview of processed samples, running jobs, and a multi-sample comparison table.
* **Samples Management:** Add new samples, including uploading custom FASTA contigs, and manage existing ones.
* **Run Pipeline:** Configure and launch different pipeline modules (QC, Assembly, Binning, etc.) per sample.
* **Results Visualization:** Interactive tables, QC/Assembly tabs, clustering charts for ARG/Virus (where available), and CSV exports.

### Starting the GUI

You can launch the GUI using the provided script:

```bash
./memamore-gui.sh
```

This will start the local Flask development server (by default on `http://127.0.0.1:5000`) and attempt to open your default web browser.

Alternatively, you can start it via the main CLI:

```bash
bin/main.sh gui --port 5000
```

---

## **Environment Management**

Each environment is created by its own script under `/envs/`, for example:

```bash
bash envs/create_env_plasmids.sh
```

Scripts handle:

* Micromamba activation
* Version-pinned installation
* Optional DB download (e.g. HOTSPOT, geNomad)

---

## **Database Overview**

| Database                  | Path                         | Purpose                     | Status                |
| ------------------------- | ---------------------------- | --------------------------- | --------------------- |
| **GTDB-Tk r226**          | `refdata/gtdbtk/release226/` | MAG taxonomy                | ✓                     |
| **PLASMe DB**             | `refdata/plasme/DB/`         | Plasmid detection           | ✓                     |
| **geNomad DB**            | `refdata/genomad/`           | Virus/plasmid detection     | ✓                     |
| **ViralVerify DB**        | `refdata/viralverify/`       | Auxiliary plasmid data      | ✓                     |
| **HOTSPOT DB**            | `external/HOTSPOT/database/` | Plasmid host prediction     | ✓                     |
| **VirSorter2 / VIBRANT / CheckV** | `refdata/{virsorter2,vibrant,checkv}/` | Virus detection and QC | ✓ (auto-download) |
| **iPHoP** (optional)      | `refdata/iphop/*_rw/`        | Virus host prediction       | ~340 GB, manual download |
| **DeepARG / CARD / MEGARes** | `refdata/{deeparg,card,amrplusplus}/` | ARG detection     | ✓ (auto-download)     |
| **VFDB / ResFinder / PlasmidFinder** | bundled with ABRicate | Virulence, ARGs, replicons | ✓                 |
| **CheckM2**               | `refdata/checkm2/`           | MAG quality                 | ✓ (auto-download)     |

---

## **Output Example**

```
SRA/
├── reads/<sample>/                     # normalized read links + manifest.tsv
├── qc/{fastp,fastqc,multiqc}/<sample>/
├── assemblies/{spades,contig_qc,metaplasmidspades,metaviralspades}/<sample>/
├── binning/{metabat2,binny,comebin,magscot,checkm2,gtdbtk,barrnap,trnascan}/<sample>/
│                                       # refined MAGs: binning/magscot/<sample>/bins/*.fa
├── plasmids/<sample>/{viralverify,genomad,plasme,mobrecon,union,cluster,typer,hotspot}/
├── viruses/<sample>/{viralverify,genomad,virsorter2,vibrant,union,checkv,iphop}/
├── args/<sample>/{deeparg,rgi,amrplusplus,abricate,summary}/
└── summary/                            # samples_overview.tsv, mags.tsv, viruses.tsv, *_matrix.tsv
```

Each module writes outputs in its own subfolder with `.done` markers to support resumable execution.

---

## **Testing (Synthetic Data)**

Two simulated test sets from public NCBI genomes (`scripts/make_mini_testdata.py`):

* `--profile mini` (default): phiX174, lambda, plasmid F, ~2 MB of reads, minutes.
  Covers ingest, QC, assembly and geNomad.
* `--profile mock`: two samples (`mockA`, `mockB`) with different abundances of
  *S. aureus*, *H. pylori*, *M. pneumoniae*, the ARG plasmid pEK499 (CTX-M-15,
  TEM-1, OXA-1, …) and phages T4/lambda/phiX; ~70 MB of reads per sample.
  Covers every module including binning and the multi-sample summary.

```bash
python3 scripts/make_mini_testdata.py --profile mock      # -> test/mock/
bin/main.sh run ingest --fastq test/mock/manifest.tsv
for s in mockA mockB; do
  for m in qc assembly binning plasmids viruses args; do bin/main.sh run $m --sample $s; done
done
bin/main.sh summary
```

Reference run (laptop, 14 cores, 15 GB RAM; ~45 min per sample, VirSorter2 ~20 min of it):
3 MAGs per sample (*H. pylori* ~100 %, *M. pneumoniae* 96–98 %, *S. aureus* 76–100 %
completeness by CheckM2), the pEK499 genes (CTX-M-15, OXA-1, TEM, sul1, tet(A), …) via
ResFinder/RGI/DeepARG, IncFII/IncFIA replicons, *H. pylori* cag/vac and *S. aureus*
virulence genes via VFDB, and lambda/phiX/T4 plus *S. aureus* prophages as viruses.
GTDB-Tk is skipped below 64 GB RAM; COMEBin fails on samples this small (too few
contigs to train), which the pipeline tolerates.

### Stub run (no environments, no databases)

`tests/stub/run.sh` runs the whole chain (ingest → qc → assembly → binning →
plasmids → viruses → args) in a few seconds with fake tools that write
correctly named outputs. It exercises the module scripts, the `.done`
bookkeeping and the GUI, not the tools themselves. It works in a throwaway
directory (`STUB_WORK`, default: a new temp dir) and never touches `SRA/` or
`refdata/`. Needs `python3` with PyYAML and `gawk`; the GUI check needs Flask.

```bash
STUB_WORK=/tmp/memamo_stub tests/stub/run.sh            # or: tests/stub/run.sh qc assembly
python3 tests/stub/check_gui.py /tmp/memamo_stub        # GUI API shows results per module
```

CI (`.github/workflows/ci.yml`) runs this on every push/PR.

---

## **Development Notes**

* Strict Bash safety (`set -euo pipefail`)
* Modular orchestration by module
* Automatic environment setup
* CUDA support where available
* CI: `bash -n`, shellcheck, flake8, GUI smoke test and the stub pipeline on every push/PR

---

## **Tool References**

A complete reference checklist for all tools used in MeMaMoRe is provided for citation in the *Methods → References* section.

### **Read Quality Control**

* Andrews S. 2010. *FastQC*.
* Ewels P *et al.* 2016. *MultiQC.* *Bioinformatics* 32(19):3047. doi:10.1093/bioinformatics/btw354.
* Chen S *et al.* 2018. *fastp.* *Bioinformatics* 34(17):i884–i890. doi:10.1093/bioinformatics/bty560.

### **Assembly**

* Prjibelski A *et al.* 2020. *SPAdes.* *Curr Protoc Bioinformatics* 70:e102.
* Nurk S *et al.* 2017. *metaSPAdes.* *Genome Res.* 27(5):824–834.
* Antipov D *et al.* 2019. *metaPlasmidSPAdes.* *Genome Res.* 29(6):961–968.
* Antipov D *et al.* 2020. *ViralVerify.* *Bioinformatics* 36(14):4126–4129.

### **Binning and Refinement**

* Kang DD *et al.* 2019. *MetaBAT2.* *PeerJ* 7:e7359.
* Hickl O *et al.* 2022. *Binny.* *Brief Bioinform.* 23(6):bbac431.
* Wang Z *et al.* 2024. *COMEBin.* *Nat Commun* 15:585.
* Rühlemann M *et al.* 2022. *MAGScoT.* *Bioinformatics* 38(24):5430–5433.

### **MAG Quality and Taxonomy**

* Chklovski A *et al.* 2023. *CheckM2.* *Nat Methods* 20:1203–1212.
* Chaumeil P-A *et al.* 2022. *GTDB-Tk.* *Bioinformatics.*
* Parks DH *et al.* 2021. *GTDB.* *Nucleic Acids Res.* 50:D785–D794.
* Bushnell B. 2014. *BBTools.* LBNL.

### **Plasmid Detection and Host Prediction**

* Tang K *et al.* 2023. *PLASMe.* *Nucleic Acids Res.* 51(15):e83.
* Robertson J, Nash JHE. 2018. *MOB-suite.* *Microb Genom* 4(8):e000206.
* Camargo AP *et al.* 2024. *geNomad.* *Nat Biotechnol* 42:1303–1312.
* Ji W *et al.* 2023. *HOTSPOT.* *Bioinformatics* 39(5):btad283.

### **Viral Analysis**

* Guo J *et al.* 2021. *VirSorter2.* *Microbiome* 9:37.
* Kieft K *et al.* 2020. *VIBRANT.* *Microbiome* 8:90.
* Nayfach S *et al.* 2021. *CheckV.* *Nat Biotechnol* 39:578–585.
* Roux S *et al.* 2023. *iPHoP.* *PLOS Biol.* 21(4):e3002083.

### **ARGs and Virulence**

* Arango-Argoty G *et al.* 2018. *DeepARG.* *Microbiome* 6:23.
* Alcock BP *et al.* 2023. *CARD.* *Nucleic Acids Res.* 51(D1):D690–D699.
* Bonin N *et al.* 2023. *MEGARes and AMR++ v3.0.* *Nucleic Acids Res.* 51(D1):D744–D752.
* Seemann T. *ABRicate.* github.com/tseemann/abricate.
* Liu B *et al.* 2022. *VFDB 2022.* *Nucleic Acids Res.* 50(D1):D912–D917.
* Florensa AF *et al.* 2022. *ResFinder.* *Microb Genom* 8(1):000748.
* Carattoli A *et al.* 2014. *PlasmidFinder.* *Antimicrob Agents Chemother* 58(7):3895–3903.

---

## **Citation**

> Sperber, S. (2025). *MeMaMoRe: A modular pipeline for metagenomic MAG, plasmid, and virus reconstruction.*
> University of Potsdam / ATB Potsdam.
