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

