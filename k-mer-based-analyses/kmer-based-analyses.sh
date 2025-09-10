#!/bin/bash

# Load conda module (adjust if needed)
module load miniforge/24.7.1-0

# Define paths (update these)
WORKFLOW_DIR="/path/to/your/snakemake/project"
SNAKEFILE="$WORKFLOW_DIR/Snakefile"
CONFIGFILE="$WORKFLOW_DIR/config/config.yaml"
ENV_PATH="/path/to/your/conda/env"  # e.g., /home/user/.conda/envs/snakemake

# Change to workflow directory
cd "$WORKFLOW_DIR" || { echo "Failed to enter workflow directory"; exit 1; }

# Activate conda environment
conda activate "$ENV_PATH"

# Run Snakemake
snakemake --cores 20 \
          --printshellcmds \
          --verbose \
          --snakefile "$SNAKEFILE" \
          --configfile "$CONFIGFILE" \
          --latency-wait 300
          # --unlock  # Uncomment if needed

# Check if Snakemake succeeded
if [ $? -ne 0 ]; then
    echo "❌ Snakemake run failed."
    conda deactivate
    exit 1
fi

# Deactivate environment
conda deactivate

echo "✅ Snakemake run completed successfully."
