#!/bin/bash

# Load required modules (make sure these are available in your environment)
module load plink/1.90beta6.16
module load bcftools/1.18

# -------- INPUT VARIABLES --------
SPECIES="SPECIES_NAME"  # Replace with your species name, e.g. RhinePhoxinus

# Prefix path to PLINK binary files (.bed, .bim, .fam)
# These files contain genotype data for your samples
BED_FILE="/path/to/plink_binary/${SPECIES}/${SPECIES}"

# Define genomic region to analyze
CHR=CHROMOSOME_NUMBER          # e.g. 3
START_POSITION=START_BP        # e.g. 47140000
END_POSITION=END_BP            # e.g. 47160000

# Define output directory and filename prefix
OUTPUT_DIR="/path/to/output/${SPECIES}/PCA/Larger_regions"
mkdir -p "${OUTPUT_DIR}"

OUTPUT_PREFIX="${OUTPUT_DIR}/PCA_kmer_suspect_${CHR}_${START_POSITION}_${END_POSITION}"

# -------- RUN PLINK --------
plink --cow --allow-no-sex \
      --nonfounders --set-missing-var-ids @:# \
      --chr ${CHR} \
      --from-bp ${START_POSITION} --to-bp ${END_POSITION} \
      --bfile ${BED_FILE} \
      --distance-matrix \
      --out ${OUTPUT_PREFIX}

# -----------------------------------------
# EXPLANATION:
# - ${BED_FILE}.bed/.bim/.fam: PLINK binary genotype files for your dataset.
# - --chr, --from-bp, --to-bp: Select SNPs from specified chromosome and base pair range.
# - --distance-matrix: Compute genetic distance matrix between samples.
# - Output:
#    * ${OUTPUT_PREFIX}.mdist - distance matrix file (for PCA/MDS)
#    * ${OUTPUT_PREFIX}.log   - PLINK log file with details of the run
# -----------------------------------------
