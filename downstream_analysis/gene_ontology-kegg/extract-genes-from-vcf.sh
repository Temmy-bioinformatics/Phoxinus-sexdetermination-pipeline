#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -M T.Oriowo@leibniz-lib.de
#$ -m n
#$ -q fast.q
#$ -e /home/toriowo/SNAKEMAKE/SexDetermination_snakemake/logs
#$ -o /home/toriowo/SNAKEMAKE/SexDetermination_snakemake/logs
#$ -pe smp 2
#$ -N Annotation

module load bedtools/2.29.2

# -------------------------
# Settings
# -------------------------
SPECIES="DanubeCsikii"
GFF_FILE="/home/toriowo/SNAKEMAKE/SexDetermination_snakemake/data/annotation/Hap1_ncbi_annotation.gff3"
REGION_FILE="/home/toriowo/SNAKEMAKE/SexDetermination_snakemake/haplome_one_results/Candidate_regions/${SPECIES}/${SPECIES}_candidate_regions.txt"
OUTPUT_DIR="/home/toriowo/SNAKEMAKE/SexDetermination_snakemake/haplome_one_results/Candidate_regions/${SPECIES}/candidate_gene_annotations_multi"

mkdir -p "$OUTPUT_DIR"

# -------------------------
# Preprocess GFF: extract gene features once into a BED-like format
# -------------------------
GENE_BED=$(mktemp)
awk -F'\t' '
  $3 == "gene" {
      chr = $1; gsub(/^Chr/, "", chr);
      split($9, info, ";");
      gene_id = "NA"; gene_name = "NA";
      for (i in info) {
          if (info[i] ~ /^ID=/) gene_id = substr(info[i], 4);
          if (info[i] ~ /^Name=/) gene_name = substr(info[i], 6);
      }
      print chr, $4-1, $5, gene_id, gene_name;
  }' OFS='\t' "$GFF_FILE" > "$GENE_BED"

# -------------------------
# Process each flanking distance
# -------------------------
for FLANK in 250000 500000 1000000; do
    echo ">>> Processing ${FLANK}bp flanking regions..."

    while read -r chr start end; do
        region_name="${chr}_${start}_${end}"
        echo "    - Region: $region_name"

        OUTPUT_FILE="${OUTPUT_DIR}/flanking_genes_${region_name}_${FLANK}.tsv"

        # Ensure flanking start coordinate >= 0
        adj_start=$((start - FLANK))
        if (( adj_start < 0 )); then adj_start=0; fi
        adj_end=$((end + FLANK))

        # Create a temporary BED file for this region
        REGION_BED=$(mktemp)
        echo -e "${chr}\t${adj_start}\t${adj_end}\t${region_name}" > "$REGION_BED"

        # Intersect genes with the region
        bedtools intersect -wa -wb -a "$GENE_BED" -b "$REGION_BED" |
            awk -F'\t' -v OFS='\t' '{print $1, $2, $3, $4, $5, $8}' > "$OUTPUT_FILE"

        rm -f "$REGION_BED"
        echo "      → Output: $OUTPUT_FILE"
    done < "$REGION_FILE"
done

# Cleanup
rm -f "$GENE_BED"
echo "All flanking regions processed successfully."
