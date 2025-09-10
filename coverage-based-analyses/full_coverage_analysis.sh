#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -M T.Oriowo@leibniz-lib.de
#$ -m n
#$ -q medium.q,large.q
#$ -e /home/toriowo/Phoxinus_SexDet_Project/Scripts/logs
#$ -o /home/toriowo/Phoxinus_SexDet_Project/Scripts/logs
#$ -pe smp 10             
#$ -t 1-6            
#$ -N Check_stats

# Load samtools and multiqc modules
module load samtools/1.19.2
module load multiqc/1.12

projectdir="/home/toriowo/Projects/SexDetermination_Project_Clean"

# Define an array of species
SPECIES_LIST=("DanubeCsikii" "RhineCsikii" "RhinePhoxinus" "SeinePhoxinus" "RhoneSeptimaniae" "WeserMorella")

# Determine the species to process based on the task ID
SPECIES=${SPECIES_LIST[$((SGE_TASK_ID - 1))]}

echo "Processing species: $SPECIES"

# Create species-specific output directories
mkdir -p ${projectdir}/Results/samtoolsstats/New/$SPECIES

# Make output directories if they don't exist yet
mkdir -p ${projectdir}/My_Scripts/logs/samtoolsstats
mkdir -p ${projectdir}/Results/samtoolsstats

# Function to print elapsed time
print_time() {
    local START=$1
    local END=$(date +%s)
    local ELAPSED=$((END - START))
    echo "Elapsed time: ${ELAPSED} seconds."
}

# Read BAM file paths from the metadata files and store them in an array
BAM_LIST=()
while IFS= read -r line; do
    BAM_LIST+=("$line")
done < <(cat ${projectdir}/Data/Sample_lists/Bamlinks/${SPECIES}_Female_Bamfiles_Metadata.txt  ${projectdir}/Data/Sample_lists/Bamlinks/${SPECIES}_Male_Bamfiles_Metadata.txt)

# Iterate over each BAM file
for bam in "${BAM_LIST[@]}"; do
    echo
    bamfile=$(basename "$bam")
    echo "Processing BAM file: $bamfile"
    echo

    out="${projectdir}/My_Scripts/logs/samtoolsstats/samtoolsstats_${bamfile}.out"
    err="${projectdir}/My_Scripts/logs/samtoolsstats/samtoolsstats_${bamfile}.err"

    # Run samtools stats
    START=$(date +%s)
    echo "Running samtools stats for file ${bam}..."
    samtools stats "$bam" > "${projectdir}/Results/samtoolsstats/New/$SPECIES/${bamfile}.samtoolsstats.txt" 2> "$err"
    print_time $START
done

# Run multiqc after all samtools processing is done
echo "Running MultiQC for species: $SPECIES"
multiqc_output_dir="${projectdir}/Results/samtoolsstats/Multi_QC_results"
mkdir -p "${multiqc_output_dir}"

# Change to MultiQC output directory
cd "${multiqc_output_dir}"

# Run multiqc on the results of the current species
multiqc "${projectdir}/Results/samtoolsstats/$SPECIES" -n "${SPECIES}_multiqc_samtools_stats"

echo "MultiQC completed for species: $SPECIES"

#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -M T.Oriowo@leibniz-lib.de
#$ -m n
#$ -q small.q,medium.q,large.q
#$ -e /home/toriowo/Projects/SexDetermination_Project_Clean/My_Scripts/logs
#$ -o /home/toriowo/Projects/SexDetermination_Project_Clean/My_Scripts/logs
#$ -t 1
#$ -pe smp 15
#$ -N Merge_1

# Load required module
module load samtools/1.19.2 || { echo "Error: Failed to load samtools module"; exit 1; }

# Define directories
output_dir="/share/pool/toriowo/Resources/Sex_Determination/Merged_Bams/Haplotype_1"
file_input_dir="/home/toriowo/Projects/SexDetermination_Project_Clean/Data/Sample_lists/Bamlinks/Updated_lists"
stats_dir="${output_dir}/samtoolsstats"

# Define arrays for species and sex
SPECIES_LIST=("RhinePhoxinus")
SEX_LIST=("Female" "Male")

# Determine the species to process based on the task ID
SPECIES=${SPECIES_LIST[$((SGE_TASK_ID - 1))]}

echo "Processing species: $SPECIES"

# Create output directories if they do not exist
mkdir -p "${output_dir}/${SPECIES}"
mkdir -p "${stats_dir}/${SPECIES}"

# Process each sex
for SEX in "${SEX_LIST[@]}"; do
    metadata_file="${file_input_dir}/${SPECIES}_${SEX}_Bamfiles_Metadata.txt"
    merged_bam="${output_dir}/${SPECIES}/${SPECIES}_${SEX}_merged.bam"
    output_file="${stats_dir}/${SPECIES}/${SPECIES}_${SEX}_merged_samtoolsstats.txt"

    # Check if metadata file exists
    if [[ ! -f "${metadata_file}" ]]; then
        echo "Warning: Metadata file not found for ${SPECIES} ${SEX}: ${metadata_file}"
        continue
    fi

    # Merge BAM files
    echo "Merging BAM files for ${SPECIES} ${SEX}..."
    samtools merge "${merged_bam}" -@ 40 -b "${metadata_file}" || { echo "Error: Failed to merge BAM files for ${SPECIES} ${SEX}"; continue; }

    # Index the merged BAM file
    echo "Indexing merged BAM file for ${SPECIES} ${SEX}..."
    samtools index "${merged_bam}" || { echo "Error: Failed to index BAM file for ${SPECIES} ${SEX}"; continue; }

    # Run samtools stats
    echo "Running samtools stats for ${SPECIES} ${SEX}..."
    samtools stats "${merged_bam}" -@ 40 > "${output_file}" || { echo "Error: Failed to run samtools stats for ${SPECIES} ${SEX}"; continue; }
done

echo "Processing complete for species: $SPECIES"

#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -M T.Oriowo@leibniz-lib.de
#$ -m n
#$ -q small.q,medium.q,large.q
#$ -e /home/toriowo/Sex_Determination/My_Scripts/logs
#$ -o /home/toriowo/Sex_Determination/My_Scripts/logs
#$ -t 1-3 
#$ -pe smp 10
#$ -N Hap2_covstats

# Load necessary module
module load samtools/1.10

# Set directories
MERGE_DIR="/share/pool/toriowo/Resources/Sex_Determination/Merged_Bams/Haplotype_2"
STATS_DIR="${MERGE_DIR}/samtoolsstats"

# List of species
SPECIES_LIST=("DanubeCsikii" "RhinePhoxinus" "WeserMorella")

# Determine the species to process based on the task ID
SPECIES=${SPECIES_LIST[$((SGE_TASK_ID - 1))]}

echo "Processing species: $SPECIES"

# Define paths to merged BAM files
FEMALE_BAM="${MERGE_DIR}/${SPECIES}/${SPECIES}_Female_merged.bam"
MALE_BAM="${MERGE_DIR}/${SPECIES}/${SPECIES}_Male_merged.bam"

# Find minimum and maximum coverage for Female reads
MIN_F=$(samtools depth ${FEMALE_BAM} | awk 'BEGIN {min=10000} {if ($3<min) min=$3} END {print min}')
MAX_F=$(samtools depth ${FEMALE_BAM} | awk 'BEGIN {max=0} {if ($3>max) max=$3} END {print max}')

# Find minimum and maximum coverage for Male reads
MIN_M=$(samtools depth ${MALE_BAM} | awk 'BEGIN {min=10000} {if ($3<min) min=$3} END {print min}')
MAX_M=$(samtools depth ${MALE_BAM} | awk 'BEGIN {max=0} {if ($3>max) max=$3} END {print max}')

# Define paths to samtools stats files
FEMALE_STATS="${STATS_DIR}/${SPECIES}/${SPECIES}_Female_merged_samtoolsstats.txt"
MALE_STATS="${STATS_DIR}/${SPECIES}/${SPECIES}_Male_merged_samtoolsstats.txt"

# Find modal coverage for Female and Male reads
MODCOV_F=$(grep ^COV ${FEMALE_STATS} | cut -f 2- | awk -v max=0 '{if ($3>max) {modal=$2; max=$3}} END {print modal}')
MODCOV_M=$(grep ^COV ${MALE_STATS} | cut -f 2- | awk -v max=0 '{if ($3>max) {modal=$2; max=$3}} END {print modal}')

# Calculate coverage ratio
COV_RATIO=$(awk -v male_cov=$MODCOV_M -v female_cov=$MODCOV_F 'BEGIN { print (male_cov / female_cov) }')

# Calculate sizes of BAM files
BAMSIZE_F=$(ls -s ${FEMALE_BAM} | awk '{print $1}')
BAMSIZE_M=$(ls -s ${MALE_BAM} | awk '{print $1}')

# Calculate BAM file size ratio
BAM_RATIO=$(awk -v male_size=$BAMSIZE_M -v female_size=$BAMSIZE_F 'BEGIN { print (male_size / female_size) }')

# Ensure the output directory for the species exists
SPECIES_OUTPUT_DIR="/home/toriowo/SNAKEMAKE/SexDetermination_snakemake/haplome_two_results/DifCover/${SPECIES}/CovStats"
mkdir -p ${SPECIES_OUTPUT_DIR}

# Initialize the output file for the species
SPECIES_OUTPUT_FILE="${SPECIES_OUTPUT_DIR}/cov_stats_${SPECIES}.tab"
echo -e "SPECIES\tf_min_cov\tf_max_cov\tm_min_cov\tm_max_cov\tcov_ratio\tbam_ratio" > ${SPECIES_OUTPUT_FILE}

# Report all results to the species-specific output file
echo -e "${SPECIES}\t${MIN_F}\t${MAX_F}\t${MIN_M}\t${MAX_M}\t${COV_RATIO}\t${BAM_RATIO}" >> ${SPECIES_OUTPUT_FILE}

echo "Processing complete for species: $SPECIES"

#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -M T.Oriowo@leibniz-lib.de
#$ -m n
#$ -q small.q,medium.q,large.q
#$ -e /home/toriowo/SNAKEMAKE/SexDetermination_snakemake/logs
#$ -o /home/toriowo/SNAKEMAKE/SexDetermination_snakemake/logs
#$ -t 1-3 
#$ -pe smp 5
#$ -N s2_DifCover

# Load necessary modules
module load bedtools/2.29.2
module load samtools/1.10
module load miniforge/24.7.1-0

conda activate /share/pool/CompGenomVert/Applications/conda-env-deeptools

# Function to print elapsed time in human-readable format
print_time() {
    local END=$(date +%s)
    local DIFF=$(($END - $1))

    local dd=$(($DIFF / 86400))
    local dh=$((($DIFF % 86400) / 3600))
    local dm=$((($DIFF % 3600) / 60))
    local ds=$(($DIFF % 60))

    if [ $dd -gt 0 ]; then
        echo "done in ${dd} days and ${dh} hours."
    elif [ $dh -gt 0 ]; then
        echo "done in ${dh} hours and ${dm} minutes."
    elif [ $dm -gt 0 ]; then
        echo "done in ${dm} minutes and ${ds} seconds."
    else
        echo "done in ${ds} seconds."
    fi
}

# List of species
SPECIES_LIST=("DanubeCsikii" "RhinePhoxinus" "WeserMorella")
SPECIES=${SPECIES_LIST[$((SGE_TASK_ID - 1))]}

# Define directories and paths
PROJECT_DIR="/home/toriowo/SNAKEMAKE/SexDetermination_snakemake"
BAM_DIR="/share/pool/toriowo/Resources/Sex_Determination/Merged_Bams/Haplotype_1/${SPECIES}"
COV_TAB="${PROJECT_DIR}/haplome_one_results/DifCover/${SPECIES}/CovStats/cov_stats_${SPECIES}.tab"
SCRIPTS_DIR="/home/toriowo/Projects/SexDetermination_Project_Clean/My_Scripts/difcover_scripts"
LOG_DIR="/home/toriowo/SNAKEMAKE/SexDetermination_snakemake/logs"
OUT_DIR="${PROJECT_DIR}/haplome_one_results/DifCover/${SPECIES}/DifCover_result"

mkdir -p "${LOG_DIR}/${SPECIES}" "${OUT_DIR}"

# Define input BAM files
BAM1="${BAM_DIR}/${SPECIES}_Female_merged.bam"
BAM2="${BAM_DIR}/${SPECIES}_Male_merged.bam"

# Extract coverage information from table
a=$(grep "${SPECIES}" "${COV_TAB}" | cut -f 2)
A=$(grep "${SPECIES}" "${COV_TAB}" | cut -f 3)
b=$(grep "${SPECIES}" "${COV_TAB}" | cut -f 4)
B=$(grep "${SPECIES}" "${COV_TAB}" | cut -f 5)
AC=$(grep "${SPECIES}" "${COV_TAB}" | cut -f 7)

v=1000  # Target number of valid bases
l=500   # Minimum window size
bin=1   # Histogram bin size
p=2     # Enrichment score threshold

# Stage 1: Convert BAMs to Union BED
cd "${OUT_DIR}"
START=$(date +%s)

# echo "# Stage 1"
# "${SCRIPTS_DIR}/from_bams_to_unionbed.sh" "${BAM1}" "${BAM2}"

# if [ ! -s "sample1_sample2.unionbedcv" ]; then
#     echo "File sample1_sample2.unionbedcv is empty, exiting."
#     exit 1
# fi
# print_time $START
# chmod 777 sample1_sample2.unionbedcv

# Stage 2: Compute ratio per window
echo "# Stage 2"
START=$(date +%s)
"${SCRIPTS_DIR}/from_unionbed_to_ratio_per_window_CC0" -a "$a" -A "$A" -b "$b" -B "$B" -v "$v" -l "$l" sample1_sample2.unionbedcv

if [ ! -s "sample1_sample2.ratio_per_w_CC0_a${a}_A${A}_b${b}_B${B}_v${v}_l${l}" ]; then
    echo "File sample1_sample2.ratio_per_w_CC0_a${a}_A${A}_b${b}_B${B}_v${v}_l${l} is empty, exiting."
    exit 1
fi

print_time $START

# Stage 3: Prepare input for DNAcopy
echo "Starting Stage 3: Prepare input for DNAcopy"
in_file="${OUT_DIR}/sample1_sample2.ratio_per_w_CC0_a${a}_A${A}_b${b}_B${B}_v${v}_l${l}"
output_file="${in_file}.log2adj_${AC}"

if [ ! -f "$in_file" ]; then
    echo "File $in_file not found, exiting."
    exit 1
fi

echo -e "scaffold\twindow_start\t${AC}*log2(ratio)" > "$output_file"
chmod u+w "$output_file"
awk -v AC="$AC" '{if($2!="window_start") print $1"\t"$2"\t"log(AC*$7)/log(2)}' "$in_file" >> "$output_file"

if [ ! -s "$output_file" ]; then
    echo "File $output_file is empty or not created, exiting."
    exit 1
fi

echo "Generated file: $output_file"
echo "Ready to run DNAcopy"
echo "Input $output_file into Rscript ${SCRIPTS_DIR}/run_DNAcopy_from_bash.R"

# Uncomment to run the R script for DNAcopy analysis
#Rscript "${SCRIPTS_DIR}/run_DNAcopy_from_bash.R" "$output_file"

# Uncomment to run the final optional R script
input_for_r="sample1_sample2.ratio_per_w_CC0_a0_A8036_b0_B8041_v1000_l500.log2adj_1.06953"
Rscript "${SCRIPTS_DIR}/run_DNAcopy_from_bash.R" "$input_for_r"
