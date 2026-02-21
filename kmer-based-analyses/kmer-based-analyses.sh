#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -M T.Oriowo@leibniz-lib.de
#$ -e /home/toriowo/PROJECTS/sex-determination-pipeline/code/job-scripts/k-mer-based-analyses/workflow/logs
#$ -o /home/toriowo/PROJECTS/sex-determination-pipeline/code/job-scripts/k-mer-based-analyses/workflow/logs
#$ -q medium.q,large.q
#$ -pe smp 21
#$ -m n
#$ -N kmer-based-sexfindr

# Print job info
echo "========================================="
echo "Job started on: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $JOB_ID"
echo "========================================="

# Load conda module (adjust if needed)
module load miniforge/24.7.1-0

# Define paths (update these)
WORKFLOW_DIR="/home/toriowo/PROJECTS/sex-determination-pipeline/code/job-scripts/k-mer-based-analyses"
SNAKEFILE="$WORKFLOW_DIR/Snakefile"
CONFIGFILE="$WORKFLOW_DIR/config/config.yaml"
ENV_PATH="/home/toriowo/.conda/envs/snakemake"

# Change to workflow directory
cd "$WORKFLOW_DIR" || { echo "Failed to enter workflow directory"; exit 1; }

# Activate conda environment
conda activate "$ENV_PATH" || { echo "Failed to activate conda environment"; exit 1; }

# Verify snakemake is available
which snakemake || { echo "Snakemake not found in environment"; exit 1; }

# Unlock if needed (uncomment if you have lock issues)
# snakemake --unlock

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
    echo "Check logs in: $WORKFLOW_DIR/workflow/logs/"
    conda deactivate
    exit 1
fi

# Deactivate environment
conda deactivate

echo "========================================="
echo "✅ Snakemake run completed successfully."
echo "Job finished on: $(date)"
echo "========================================="