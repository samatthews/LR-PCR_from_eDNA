#!/bin/bash
# Step 0 - install the tools into a conda environment called "lrpcr".
# Needs conda (Miniforge recommended: https://github.com/conda-forge/miniforge).
# Run once:   bash 00_setup_env.sh      then:   conda activate lrpcr
set -euo pipefail

PKGS="blast>=2.14 entrez-direct seqkit taxonkit python=3.11 pandas numpy"

if conda env list | grep -q "^lrpcr "; then
  echo "Environment 'lrpcr' already exists. To rebuild: conda env remove -n lrpcr"
  exit 0
fi

# Try a native build first. On Apple Silicon some bioconda packages are Intel-only;
# if the native solve fails, fall back to Intel builds (they run under Rosetta).
if ! conda create -y -n lrpcr -c conda-forge -c bioconda $PKGS; then
  echo "Native install failed; retrying with Intel (osx-64) packages..."
  CONDA_SUBDIR=osx-64 conda create -y -n lrpcr -c conda-forge -c bioconda $PKGS
  conda run -n lrpcr conda config --env --set subdir osx-64
fi

echo
echo "Done. Check with:"
echo "  conda activate lrpcr && blastn -version && efetch -version && seqkit version && taxonkit version"
