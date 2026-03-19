#!/usr/bin/env bash
set -euo pipefail
micromamba create -n amrplusplus_env -c bioconda -c conda-forge bwa samtools python bedtools -y
