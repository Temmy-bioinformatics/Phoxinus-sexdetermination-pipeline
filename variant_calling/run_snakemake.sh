#!/bin/bash
set -euo pipefail

# Change to your project directory (replace with your path)
cd /path/to/your-project

# Activate conda base or ensure conda is in PATH
# source ~/miniconda3/etc/profile.d/conda.sh

# (Optional) Activate your snakemake conda environment if you want
# conda activate snakemake

# Run Snakemake with conda environments and 8 CPU cores (adjust cores as needed)
snakemake --cores 8 --use-conda --printshellcmds --verbose --reason --latency-wait 60

# (Optional) Deactivate conda environment if you activated it above
# conda deactivate
