#!/usr/bin/env bash
set -euo pipefail
micromamba create -n rgi_env -c bioconda -c conda-forge rgi -y
