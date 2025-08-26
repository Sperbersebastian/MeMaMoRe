````markdown
# MeMaMoRe – Metagenomic Modular Research Pipeline

This project provides a **modular pipeline** for metagenomic analyses.  
Each module is isolated, runs inside its own Micromamba environment, and can be combined or executed individually.

---

## Installation

1. Clone the repo:
   ```bash
   git clone https://github.com/Sperbersebastian/MeMaMoRe.git
   cd MeMaMoRe
````

2. Ensure [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) is installed and on your `PATH`.

3. Create the environment for the **ingest** module (other modules will add their own envs):

   ```bash
   micromamba create -y -n env_sra_tools -c conda-forge -c bioconda python=3.10 pyyaml
   ```

---

## Usage

All execution is routed through the **main driver**:

```bash
bin/main.sh <command> <module> [options]
```

### Commands

* `run` – execute a given module
* `help` – print usage info

### Ingest module (FASTQ)

Prepare raw FASTQs (from SRA download or local files) and optional 10% testsets:

```bash
bin/main.sh run ingest --fastq config/fastq_manifest.tsv
```

* Input: `config/fastq_manifest.tsv`
  Example structure:

  ```
  sample_id    country    management    R1    R2
  testA        ISR        Conv          /path/to/testA_1.fastq.gz   /path/to/testA_2.fastq.gz
  testB        USA        Org           /path/to/testB_1.fastq.gz   /path/to/testB_2.fastq.gz
  ```

  * `R2` may be left empty for single-end samples.
  * If the manifest is malformed, the pipeline will fail early. Ensure all columns exist.

* Output:

  * Raw FASTQs → `SRA/raw_fastq/<sample>/`
  * Testset FASTQs (10% downsampled) → `SRA/testsets/<sample>/`
  * Sample table → `config/samples.tsv`
  * Event log → `api/status.jsonl`

### Status + Logs

Check recent module progress:

```bash
tail -n 10 api/status.jsonl | jq .
```

Inspect logs:

```bash
ls logs/
```

---

## Roadmap (Modules)

* [x] **ingest/fastq** – normalize input & testsets
* [ ] **qc** – fastp → fastqc → multiqc
* [ ] **assembly** – metaSPAdes
* [ ] **binning** – COMEBin, Binny, MetaBAT2
* [ ] **refinement** – MAGScoT
* [ ] **taxonomy** – GTDB-Tk, BBTools
* [ ] **plasmid detection** – geNomad, PlasMe, MOB-suite
* [ ] **virus detection** – VirSorter2, VIBRANT, CheckV
* [ ] **ARGs/virulence** – DeepARG, RGI, VFDB
* [ ] **diversity + stats** – R-based analyses
* [ ] **GUI layer (later)**

---

## Contributing

* Keep each module self-contained (`modules/<name>/<tool>/exec.sh`)
* Environments are **tool-specific** (`env_<tool>`)
* Always log to `logs/` and events to `api/status.jsonl`

---


```markdown
# 🍏 MeMaMoRe: Metagenomic Analysis of Mobile Genetic Elements in Apple Orchards

This repository contains a reproducible pipeline for analyzing metagenomic shotgun data from conventional and organic apple orchards.  
The pipeline identifies and characterizes **contigs, MAGs, plasmids, viruses, and ARGs**.

---

## 📂 Project Structure

```

bin/                # Main entry scripts
configs/            # Config files
modules/            # Individual pipeline modules
SRA/                # Output data (assemblies, bins, QC, etc.)
logs/               # Runtime logs

````

---

## 🔧 Dependencies

All software is installed in isolated **micromamba environments**.  
Each module automatically activates the correct environment when run.

Main tools:

- **Fastp, FastQC, MultiQC** → read QC  
- **SPAdes/metaSPAdes** → assembly  
- **BWA-MEM2** → mapping  
- **samtools** → BAM handling  
- **CoverM** → coverage statistics  
- **CheckM2, GUNC** → MAG QC  
- **geNomad, PlasMe, MOB-suite** → plasmid detection  
- **VirSorter2, VIBRANT, CheckV** → virus detection  
- **DeepARG, VFDB, RGI** → resistance/virulence factors  

---

## 🏃 Usage

Run modules via the main entry script:

```bash
bin/main.sh run <module> [--sample SAMPLE] [--test] [--force]
````

Options:

* `--test` → runs in quick test mode (subset of reads, reduced runtime)
* `--force` → overwrites existing outputs

---

## A. Quality Control

Raw reads are trimmed and quality-checked:

* **fastp** → adapter removal, trimming, filtering
* **FastQC** → per-sample QC
* **MultiQC** → combined summary

Outputs under:

```
SRA/fastp/
  ├── *_trimmed.fastq.gz
  └── fastqc_results/
```

---

## B. Assembly (SPAdes)

Reads are assembled with **SPAdes/metaSPAdes** into contigs.

Outputs under:

```
SRA/assemblies/spades/<sample>/
  ├── contigs.fasta
  ├── contigs.len1000.fasta      # contigs ≥ 1000 bp
  ├── logs/
  └── map/                       # BAMs + coverage stats
```

---

## C. Assembly QC (contig\_qc)

After assembly, contigs are filtered based on mapping support.

### 🚀 Run

```bash
# Test run (subset of reads)
bin/main.sh run assembly --test --sample SAMPLE_NAME --force

# Full run
bin/main.sh run assembly --sample SAMPLE_NAME --force
```

### 📂 Outputs

Results for each sample:

```
SRA/assemblies/
├── spades/<sample>/            # raw assembly + mapping
│   ├── contigs.len1000.fasta
│   └── map/sample.sorted.bam
│
└── contig_qc/<sample>/         # QC outputs
    ├── contigs_kept.tsv         # kept contigs
    ├── contigs_dropped.tsv      # dropped contigs + reason
    ├── contigs.filtered.fasta   # FASTA with kept contigs
    └── mapping_counts.tsv       # read support summary
```

### 📊 Filtering rules

* **Breadth ≥ 0.8** → ≥80% of bases covered
* **Mean depth ≥ 2** → average coverage ≥ 2×

Contigs failing these are dropped with reason (`depth`, `breadth`).

### 🧮 Mapping Counts (`mapping_counts.tsv`)

| sample | total\_primary | mapped\_unfiltered | mapped\_filtered | prop\_unfiltered | prop\_filtered |
| ------ | -------------- | ------------------ | ---------------- | ---------------- | -------------- |
| testA  | 4,228,836      | 3,181,216          | 3,180,648        | 0.752268         | 0.752133       |

* **total\_primary** → all primary reads in BAM (excl. secondary/supplementary)
* **mapped\_unfiltered** → reads mapped to all contigs
* **mapped\_filtered** → reads mapped to QC-kept contigs only
* **prop**\* → fraction of mapped reads

---

## D. Binning

(coming soon — COMEBin, Binny, MetaBAT2 → MAGs, refinement with MAGScoT, QC with CheckM2, GUNC)

---

## E. Plasmid Detection

(coming soon — geNomad, PlasMe, MOB-suite, viralverify+genomad merge)

---

## F. Virus Detection

(coming soon — VirSorter2, VIBRANT, CheckV, deepPHAGE)

---

## G. Functional Annotation

(coming soon — DeepARG, VFDB, RGI, CARD integration)

---

## 📝 Notes

* All modules log to `logs/`
* Test mode (`--test`) is recommended for pipeline debugging
* Full runs are compute-intensive; use SLURM or HPC batch jobs

---
D. Binning (MAG recovery)

Input:
Filtered contigs from C. contig_qc.

Steps:

Initial binning
Run three independent binners in parallel:

COMEBin (SRA/binning/comebin/<sample>/)

Binny (SRA/binning/binny/<sample>/)

MetaBAT2 (SRA/binning/metabat2/<sample>/)

Each produces draft bins as FASTA files.

Bin refinement (MAGScoT)

Combine outputs of the three binners.

Run MAGScoT to merge, dereplicate, and refine.

Output: SRA/binning/magscot/<sample>/bins/ (final MAG FASTAs).

Quality control
Each MAG is evaluated with:

CheckM2 → genome completeness & contamination.

barnap → rRNA detection.

tRNAscan-SE → tRNA detection.

Recommended filters:

High-quality MAGs: ≥90% completeness, ≤5% contamination, ≥18 tRNAs, rRNA detected.

Medium-quality MAGs: ≥50% completeness, ≤10% contamination.

All others discarded.

Taxonomy assignment

Run GTDB-Tk on refined, QC-passed MAGs.

Output taxonomy tables + annotated MAG FASTAs: SRA/binning/gtdbtk/<sample>/.

Outputs per sample:

binning/comebin/ → raw COMEBin bins

binning/binny/ → raw Binny bins

binning/metabat2/ → raw MetaBAT2 bins

binning/magscot/ → refined MAGs

binning/checkm2/ → QC metrics

binning/barnap/ + binning/trnascan/ → rRNA/tRNA results

binning/gtdbtk/ → taxonomy assignment