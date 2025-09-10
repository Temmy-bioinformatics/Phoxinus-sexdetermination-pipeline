#!/bin/bash
#
# Job script template for mitogenome phylogeny pipeline
#
# Scheduler directives (SGE example; modify for your scheduler):
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -m n
#$ -q medium.q,large.q
#$ -pe smp 15
#$ -e <PATH_TO_LOG_DIR>
#$ -o <PATH_TO_LOG_DIR>
#$ -N mitogenome_phylo
#$ -M your.email@domain.com

set -euo pipefail
IFS=$'\n\t'

# === LOAD MODULES ===
module load iqtree/2.2.2.7
module load bcftools/1.19
module load mafft/7.453we
module load samtools/1.19.2
module load seqtk/1.3

# === USER CONFIGURATION ===

# Base directory where MitoExtractor output with *_final_genes_NT.fasta files are stored
OUTPUT_BASE="<PATH_TO_MITOEXTRACTOR_OUTPUT>"

# Directory where pipeline will run and intermediate files be stored
WORKING_DIR="<PATH_TO_PIPELINE_WORKDIR>"

# Path to AMAS.py script (can be a conda env executable or system path)
AMAS_PATH="<PATH_TO_AMAS.py>"

# List of genes to extract and align
GENES_OF_INTEREST=("COX1" "CYTB" "ND2" "ND4" "ND5" "ATP6" "ATP8" "COX2" "COX3" "ND1" "ND3" "ND6")

# Number of CPU threads for IQ-TREE
THREADS=15

# Outgroup taxon name for IQ-TREE - adjust as needed
OUTGROUP="SRR18560773"

# === SETUP ===
mkdir -p "$WORKING_DIR/gene_alignments"
cd "$WORKING_DIR"

# === STEP 1: Collect all *_final_genes_NT.fasta files ===
echo "Collecting gene fasta files..."
find "$OUTPUT_BASE" -type f -name "*_final_genes_NT.fasta" > all_final_genes.txt

# === STEP 2: Extract per-gene fasta sequences with sample headers ===
echo "Extracting per-gene FASTAs with sample headers..."
rm -rf gene_alignments/*
while read -r file; do
    sample=$(basename "$file" | cut -d'_' -f1)
    for gene in "${GENES_OF_INTEREST[@]}"; do
        awk -v gene="$gene" -v sample="$sample" '
            BEGIN {found=0}
            /^>/ {
                found=0
                if ($0 ~ gene) {
                    found=1
                    print ">" sample
                }
            }
            !/^>/ {
                if (found) print
            }
        ' "$file" >> "gene_alignments/${gene}.fasta"
    done
done < all_final_genes.txt

# === STEP 3: Align each gene with MAFFT ===
echo "Aligning genes with MAFFT..."
for gene in "${GENES_OF_INTEREST[@]}"; do
    mafft --localpair --maxiterate 1000 "gene_alignments/${gene}.fasta" > "gene_alignments/${gene}_aligned.fasta"
done

# === STEP 4: Concatenate alignments into supermatrix using AMAS ===
echo "Concatenating alignments into supermatrix..."
cd gene_alignments
"${AMAS_PATH}" concat -i *_aligned.fasta -f fasta -d dna -u phylip -t supermatrix.phy -p partitions.txt

# Modify partitions file for IQ-TREE format
sed 's/^\(p[0-9]*_[^ ]*\)_aligned[ ]*= */DNA, \1 = /' partitions.txt > partitions_iqtree.txt

# === STEP 5: Run IQ-TREE phylogeny inference ===
echo "Running IQ-TREE..."
iqtree2 -s supermatrix.phy \
       -p partitions_iqtree.txt \
       -m MFP+MERGE \
       -B 1000 \
       -alrt 1000 \
       -o "$OUTGROUP" \
       -T AUTO \
       -ntmax "$THREADS"

echo "Pipeline complete. Tree file is located at:"
realpath supermatrix.phy.treefile
