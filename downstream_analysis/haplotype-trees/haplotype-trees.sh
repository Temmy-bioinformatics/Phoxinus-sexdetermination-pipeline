#!/bin/bash

# Usage: ./script.sh <SPECIES_NAME>
# Example: ./script.sh RhinePhoxinus

set -euo pipefail

SPECIES=${1:?Species name must be provided, e.g. RhinePhoxinus}

# === CONFIGURE BASE PATHS ===
BASE_DIR="/path/to/SexDetermination_snakemake/haplome_one_results"
WORK_DIR="${BASE_DIR}/${SPECIES}"

PHENO_FILE="${WORK_DIR}/SNP_Gwas/pheno_files/${SPECIES}/${SPECIES}.pheno.txt"
REGIONS_FILE="${WORK_DIR}/Candidate_regions/${SPECIES}/${SPECIES}_candidate_regions.txt"
REFERENCE_FASTA="/path/to/Final_Phoxy_1_rearranged.fasta"
BCF_FILE="/path/to/Variant_calling_pipeline/haplome_one_results/${SPECIES}/phased/${SPECIES}_phased_NoMissing.recode_fixed.vcf.gz"

VCF_DIR="${WORK_DIR}/individual_bcf/${SPECIES}_Regions/combined_vcf"
FASTA_DIR="${WORK_DIR}/individual_bcf/${SPECIES}_Regions/fasta"
TREE_DIR="${WORK_DIR}/individual_bcf/${SPECIES}_Regions/haplotype_fixed_trees"

mkdir -p "$VCF_DIR" "$FASTA_DIR" "$TREE_DIR"

# === Load necessary modules (adjust or remove if not using a module system) ===
module load vcftools/0.1.16
module load samtools/1.10
module load htslib/1.10.2
module load bcftools/1.10.2
module load iqtree/2.2.2.7 

# === Map sample IDs to gendered names ===
declare -A SAMPLE_MAP
while IFS=$'\t' read -r FID IID GENDER; do
    if [[ $GENDER -eq 1 ]]; then
        SAMPLE_MAP[$IID]="Male_${IID}"
    elif [[ $GENDER -eq 2 ]]; then
        SAMPLE_MAP[$IID]="Female_${IID}"
    fi
done < "$PHENO_FILE"

# === Process each candidate region ===
while IFS=$'\t' read -r CHROM START END; do
    REGION_VCF_DIR="$VCF_DIR/Chr${CHROM}_${START}_${END}"
    REGION_FASTA_DIR="$FASTA_DIR/Chr${CHROM}_${START}_${END}"
    REGION_TREE_DIR="$TREE_DIR/Chr${CHROM}_${START}_${END}"
    FULL_VCF_DIR="$VCF_DIR/full/Chr${CHROM}_${START}_${END}"

    mkdir -p "$REGION_VCF_DIR" "$REGION_FASTA_DIR" "$REGION_TREE_DIR" "$FULL_VCF_DIR"

    for SAMPLE in "${!SAMPLE_MAP[@]}"; do
        GENDERED_NAME="${SAMPLE_MAP[$SAMPLE]}"
        
        # Extract sample VCF for the region
        bcftools view -s "$SAMPLE" -r "Chr${CHROM}:${START}-${END}" -Oz \
            -o "${REGION_VCF_DIR}/${GENDERED_NAME}_Chr${CHROM}_${START}_${END}.vcf.gz" "$BCF_FILE"
        bcftools index "${REGION_VCF_DIR}/${GENDERED_NAME}_Chr${CHROM}_${START}_${END}.vcf.gz"

        # Create haplotype consensus fasta for hap1 and hap2
        samtools faidx "${REFERENCE_FASTA}" "Chr${CHROM}:${START}-${END}" | \
        bcftools consensus -H 1 "${REGION_VCF_DIR}/${GENDERED_NAME}_Chr${CHROM}_${START}_${END}.vcf.gz" | \
        awk -v name="$GENDERED_NAME" '{if(substr($0,1,1)==">") print ">"name"_hap1_"substr($0,2); else print $0}' \
        > "${REGION_FASTA_DIR}/${GENDERED_NAME}_hap1.fa"

        samtools faidx "${REFERENCE_FASTA}" "Chr${CHROM}:${START}-${END}" | \
        bcftools consensus -H 2 "${REGION_VCF_DIR}/${GENDERED_NAME}_Chr${CHROM}_${START}_${END}.vcf.gz" | \
        awk -v name="$GENDERED_NAME" '{if(substr($0,1,1)==">") print ">"name"_hap2_"substr($0,2); else print $0}' \
        > "${REGION_FASTA_DIR}/${GENDERED_NAME}_hap2.fa"
    done

    # Extract full region VCF across all samples
    bcftools view -r "Chr${CHROM}:${START}-${END}" -Oz \
        -o "${FULL_VCF_DIR}/Chr${CHROM}_${START}_${END}.vcf.gz" "$BCF_FILE"

    # Concatenate all haplotype fasta files for tree building
    cat "${REGION_FASTA_DIR}"/*.fa > "${REGION_FASTA_DIR}/${SPECIES}_Chr${CHROM}_${START}_${END}_all_seqs.fasta"

    # Run IQ-TREE for phylogeny on concatenated sequences
    iqtree -s "${REGION_FASTA_DIR}/${SPECIES}_Chr${CHROM}_${START}_${END}_all_seqs.fasta" \
           --mem 500G -T 8 -B 1000 \
           -pre "${REGION_TREE_DIR}/tree_${SPECIES}_Chr${CHROM}_${START}_${END}"

done < "$REGIONS_FILE"

echo "All regions and samples for ${SPECIES} have been processed successfully."
