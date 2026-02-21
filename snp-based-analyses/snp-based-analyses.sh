#!/bin/bash

# Load conda module (adjust if needed)
module load miniforge/24.7.1-0

# Define paths (update these)
WORKFLOW_DIR="../snp-based-analyses"
SNAKEFILE="$WORKFLOW_DIR/Snakefile"
CONFIGFILE="$WORKFLOW_DIR/config/config.yaml"
ENV_PATH="/home/toriowo/.conda/envs/snakemake"

# Change to workflow directory
cd "$WORKFLOW_DIR" || { echo "Failed to enter workflow directory"; exit 1; }

# Activate conda environment
conda activate "$ENV_PATH"

# Run Snakemake with --use-conda flag (CRITICAL!)
snakemake --cores 20 \
          --use-conda \
          --printshellcmds \
          --verbose \
          --snakefile "$SNAKEFILE" \
          --configfile "$CONFIGFILE" \
          --latency-wait 300 \
          --rerun-incomplete
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