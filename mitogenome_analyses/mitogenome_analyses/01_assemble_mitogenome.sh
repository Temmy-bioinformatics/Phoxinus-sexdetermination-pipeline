#!/bin/bash
#
# Job submission script template for running MitoFinder on species-sex fastq datasets with job array support
#
# Instructions:
# - Replace all placeholders (e.g., <PATH_TO_...>) with actual paths before running.
# - Adjust job scheduler directives if needed.
# - Submit with your cluster's job submission command (e.g., qsub or sbatch).
#
# SGE directives (modify or remove if using another scheduler):
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -m n
#$ -e <PATH_TO_LOG_DIR>
#$ -o <PATH_TO_LOG_DIR>
#$ -q medium.q,large.q
#$ -pe smp 10
#$ -t 1-4       # Adjust number of tasks to match your species-sex combos
#$ -N Mitofinder_loop
#$ -M your.email@domain.com

set -euo pipefail
IFS=$'\n\t'

# Load required modules
module load java/jre1.8.0_231

# Optionally source your bashrc or environment variables if needed
# source ~/.bashrc

log() {
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] $*"
}

### User configurable variables ###

# Number of threads per job
THREADS=10

# Root directory containing species and sex subdirectories with fastq files
SOURCE_DIR="<PATH_TO_FASTQ_ROOT_DIR>"

# Temporary output directory for intermediate MitoFinder results
DEST_DIR="<PATH_TO_MITOFINDER_OUTPUT_DIR>"

# Final output directory to store cleaned mitochondrial contigs
FINAL_OUTPUT_DIR="<PATH_TO_FINAL_OUTPUT_DIR>"

# Reference mitochondrial genome file (GenBank or FASTA format)
REFERENCE="<PATH_TO_REFERENCE_MITOCHONDRIA>"

# Define species and sex combinations corresponding to your array job indices
SPECIES=("Species1" "Species2")
SEX=("Male" "Female")

# Calculate zero-based index from SGE_TASK_ID environment variable
task_id=$((SGE_TASK_ID - 1))

# Determine species and sex indices for this task
species_index=$((task_id / 2))
sex_index=$((task_id % 2))

species=${SPECIES[$species_index]}
sex=${SEX[$sex_index]}

log "Starting task $SGE_TASK_ID: Processing species='$species', sex='$sex'"

# Path for species-sex specific fastq files
species_sex_dir="${SOURCE_DIR}/${species}/${sex}"

# Find all *_1.fq.gz or *_1.fastq.gz files for this species-sex combo
fastq_files=($(find "$species_sex_dir" -type f \( -name "*_1.fq.gz" -o -name "*_1.fastq.gz" \)))

for fq1 in "${fastq_files[@]}"; do
    # Derive matching _2 file path by replacing _1 with _2
    fq2=$(echo "$fq1" | sed -E 's/_1\.(fq|fastq)\.gz$/_2.\1.gz/')

    # Check if matching mate file exists
    if [[ ! -f "$fq2" ]]; then
        log "Warning: Mate pair file not found for $fq1, skipping sample."
        continue
    fi

    # Extract sample name without read pair suffix and extensions
    sample=$(basename "$fq1" | sed -E 's/_(1|2)\.(fq|fastq)\.gz//')

    # Define output directory for this sample's intermediate results
    sample_out_dir="${DEST_DIR}/${species}/${sex}/${sample}"

    if [[ -d "$sample_out_dir" && -n $(find "$sample_out_dir" -mindepth 1 -type f) ]]; then
        log "Skipping $sample: Output directory already exists and contains files."
        continue
    fi

    mkdir -p "$sample_out_dir"
    cd "$sample_out_dir" || { log "ERROR: Cannot cd to $sample_out_dir"; exit 1; }

    log "Running MitoFinder for sample $sample..."
    mitofinder -j "$sample" \
               -1 "$fq1" \
               -2 "$fq2" \
               -r "$REFERENCE" \
               -o "$sample_out_dir" \
               -p "$THREADS" \
               -o 2

    # Path to expected final mtDNA contig output file
    final_mito_contig="${sample_out_dir}/${sample}_MitoFinder_megahit_mitfi_Final_Results/${sample}_mtDNA_contig.fasta"

    if [[ -f "$final_mito_contig" ]]; then
        # Destination path for cleaned output
        final_destination="${FINAL_OUTPUT_DIR}/${species}/${sex}/${sample}_mtDNA_contig.fasta"
        mkdir -p "$(dirname "$final_destination")"
        mv "$final_mito_contig" "$final_destination"
        log "Moved mtDNA contig to $final_destination"

        # Clean up intermediate output to save space
        rm -rf "$sample_out_dir"
        log "Removed intermediate files for sample $sample"
    else
        log "ERROR: Expected output file not found for sample $sample: $final_mito_contig"
    fi
done

log "Task $SGE_TASK_ID complete for species='$species', sex='$sex'."
