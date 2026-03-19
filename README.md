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
| **Plasmids**     | Detection, clustering, and host prediction           | PLASMe, MOB-suite, geNomad (✓ core), ViralVerify (✓ core), HOTSPOT (⚠ partial) |
| **Viruses** ⚠    | Viral prediction and QC *(not implemented yet)*      | geNomad, VirSorter2, VIBRANT, CheckV, iPHoP                                    |
| **Functional** ⚠ | ARG and virulence detection *(not implemented yet)*  | DeepARG, VFDB, CARD                                                            |
| **Simulation** ⚠ | Synthetic reads and controlled variation *(partial)* | ART, custom `sim_env`                                                          |

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
│       ├── viruses.sh       # ⚠ placeholder for viral analysis
│       ├── functional.sh    # ⚠ placeholder for ARG/virulence annotation
│       └── sim.sh           # Simulation orchestrator
├── envs/                    # Micromamba environment creation scripts
│   ├── create_env_plasmids.sh
│   ├── create_env_binning.sh
│   ├── create_env_qc.sh
│   └── ...
├── external/                # External cloned tools (e.g. Binny, COMEBin, HOTSPOT)
├── modules/                 # Tool-specific executors (assembly, plasmids, binning)
├── refdata/                 # Databases (GTDB, PLASMe, geNomad, HOTSPOT, etc.)
├── SRA/                     # Main workspace for input/output per sample
│   ├── reads/
│   ├── assemblies/
│   ├── binning/
│   ├── plasmids/
│   ├── viruses/
│   ├── functional/
│   └── qc/
└── config/                  # YAML configurations (params, testsets)
```

⚠ Modules **viruses**, **functional**, and **simulation** are placeholders or partially implemented.

---

## **Installation**

### 1. Prerequisites

* Linux (tested on Ubuntu 22.04 and CentOS 7)
* Micromamba or Conda
* Git ≥ 2.30
* Bash ≥ 5.0

### 2. Clone the Repository

```bash
git clone https://github.com/Sperbersebastian/MeMaMoRe.git
cd MeMaMoRe
```

### 3. Create Environments

```bash
bin/main.sh env create plasmids --force
```

Installs:

* `plasmids_core`
* `plasme_env`
* `mobsuite_env`
* `hotspot_env`
* `genomad_env`
* `viralverify_env`

⚠ `bin/main.sh env create all` not yet implemented.

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
| **HOTSPOT DB**            | `refdata/HOTSPOT/database/`  | Plasmid host prediction     | ⚠ partial             |
| **CARD / DeepARG / VFDB** | `refdata/functional/`        | ARG and virulence databases | ⚠ pending integration |

---

## **Output Example**

```
SRA/
└── test_sample/
    ├── qc/
    ├── assemblies/spades/contigs.fasta
    ├── binning/comebin/
    ├── plasmids/plasme_hp/
    ├── plasmids/mobsuite/
    ├── viruses/genomad/        # ⚠ placeholder
    └── functional/deeparg/     # ⚠ placeholder
```

Each module writes outputs in its own subfolder with `.done` markers to support resumable execution.

---

## **Testing (Synthetic Data)**

Quick test run using internal mini dataset:

```bash
bin/main.sh run sim --set miniset
bin/main.sh run plasmids --sample miniset
```

⚠ Simulation module partially implemented.

---

## **Development Notes**

* Strict Bash safety (`set -euo pipefail`)
* Modular orchestration by module
* Automatic environment setup
* CUDA support where available
* ⚠ Not yet available:

  * Virus and functional modules
  * Multi-sample summary
  * Combined “env create all”

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

### **Viral Analysis (Planned)**

* Guo J *et al.* 2021. *VirSorter2.* *Microbiome* 9:37.
* Kieft K *et al.* 2020. *VIBRANT.* *Microbiome* 8:90.
* Nayfach S *et al.* 2021. *CheckV.* *Nat Biotechnol* 39:578–585.
* Roux S *et al.* 2023. *iPHoP.* *PLOS Biol.* 21(4):e3002083.

### **Functional Annotation (Planned)**

* Arango-Argoty G *et al.* 2018. *DeepARG.* *Microbiome* 6:23.
* Alcock BP *et al.* 2023. *CARD.* *Nucleic Acids Res.* 51(D1):D690–D699.

---

## **Citation**

> Sperber, S. (2025). *MeMaMoRe: A modular pipeline for metagenomic MAG, plasmid, and virus reconstruction.*
> University of Potsdam / ATB Potsdam.
