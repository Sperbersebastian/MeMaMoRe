#!/usr/bin/env bash
set -euo pipefail
# RGI >= 6.0.5 calls "diamond --xml-blord-format", which DIAMOND 0.8.x (still allowed by
# the bioconda recipe, "diamond >=0.8.36") does not know -> pin a current DIAMOND.
micromamba create -n rgi_env -c bioconda -c conda-forge rgi "diamond>=2.1.10" -y
