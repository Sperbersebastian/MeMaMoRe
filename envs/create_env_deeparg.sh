#!/usr/bin/env bash
set -euo pipefail
micromamba create -n deeparg_env -c bioconda -c conda-forge deeparg -y
