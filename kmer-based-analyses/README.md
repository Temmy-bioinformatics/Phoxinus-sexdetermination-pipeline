K-mer Based Sex Determination Pipeline – README
Overview

This pipeline performs a k-mer based genome-wide association study (GWAS) to identify k-mers associated with sex in a given species. The workflow includes k-mer counting, strand annotation, kinship correction, GWAS association testing, extraction of significant k-mers, assembly/BLAST of top k-mers, and final visualization.

It is modular and controlled via a Snakemake master workflow, with toggleable steps for flexibility.

Project Structure
PROJECT_DIR/
├── config/
│   └── config.yaml              # Main configuration file
├── workflow/
│   ├── rules/
│   │   ├── 01_files_prep.smk
│   │   ├── 02_kmc.smk
│   │   ├── 03_strand_info.smk
│   │   ├── 04_list_kmers.smk
│   │   ├── 05_kinship_table.smk
│   │   ├── 06_kmer_gwas.smk
│   │   ├── 07_blast_kmers.smk
│   │   └── 09_plot_results.smk
│   └── scripts/
│       ├── plink_to_abyss_kmers.py
│       └── run_plot_kmer_blast.R
├── {SPECIES}/                    # Per-species working directory
│   ├── {sample}/                 # Individual samples
│   │   ├── *_1.fq.gz
│   │   ├── *_2.fq.gz
│   │   ├── output_kmc_*.kmc_pre
│   │   └── kmers_with_strand
│   ├── kmers_list_paths_final.txt
│   ├── kmers_to_use
│   ├── kmers_table.table
│   ├── kmers_table.kinship
│   ├── kmers.assoc
│   ├── kmers.assoc.tab
│   ├── most_assoc_kmers.list
│   └── blast/
│       ├── kmers_blast.out
│       └── mean_contig_len
└── plots/
    └── {SPECIES}_kmer_counts_per_chr.{png,pdf,svg}

Pipeline Modules
01 – Files Preparation (01_files_prep.smk)

Copies raw FASTQ files from the source directory to the project directory.

Handles paired-end samples automatically.

Output: {sample}_1.fq.gz and {sample}_2.fq.gz in the species directory.

02 – K-mer Counting (02_kmc.smk)

Uses KMC to count k-mers in all samples.

Generates both canonical and all k-mers counts.

Output: output_kmc_canon.kmc_pre, output_kmc_all.kmc_pre.

03 – Add Strand Information (03_strand_info.smk)

Adds strand information to k-mer counts.

Uses kmers_add_strand_information binary.

Output: kmers_with_strand.

04 – List K-mers (04_list_kmers.smk)

Combines k-mers from all samples.

Applies MAC (minimum allele count) filtering and removes low-frequency k-mers.

Outputs:

kmers_list_paths_final.txt – paths to per-sample k-mers

dirlist.txt – list of sample directories

kmers_to_use – list of filtered k-mers for GWAS

05 – Kinship Table (05_kinship_table.smk)

Builds a presence/absence k-mer table for all samples.

Computes kinship matrix to correct for population structure using emma_kinship_kmers.

Outputs:

kmers_table.table

kmers_table.kinship

06 – K-mer GWAS (06_kmer_gwas.smk)

Performs GWAS using PLINK with kinship correction.

Converts k-mer table to PLINK format using kmers_table_to_bed.

Outputs:

kmers.assoc – PLINK association results

kmers.assoc.tab – tab-delimited, human-readable table

07 – Extract Significant K-mers (Terminal Step)

Finds most significant k-mers by p-value.

Filters top k-mers for downstream analysis.

Commands:

cut -f 9 kmers.assoc.tab | grep -v 'P' | sort -g | head -1
awk '$9 <= 3.414e-07' kmers.assoc.tab > most_significant_assoc.tab
cat most_significant_assoc.tab | cut -f 2 > most_assoc_kmers.list


Output:

most_significant_assoc.tab

most_assoc_kmers.list

08 – BLAST Top K-mers (07_blast_kmers.smk)

Assembles k-mers using ABYSS.

Maps assembled sequences to reference genome with BLAST.

Outputs:

assoc_kmers_presab.txt

plink_abyss_input.txt

plink_abyss_output.txt

blast/kmers_blast.out

blast/mean_contig_len

09 – Plot Results (09_plot_results.smk)

Generates chromosome-wide visualizations of top k-mers.

Uses R script run_plot_kmer_blast.R.

Outputs:

{SPECIES}_kmer_counts_per_chr.png

{SPECIES}_kmer_counts_per_chr.pdf

{SPECIES}_kmer_counts_per_chr.svg

plots/{SPECIES}_kmer_blast_mapping_summary.csv

Configuration (config/config.yaml)

Key parameters include:

PROJECT_DIR – base working directory

SOURCE_DIRECTORY – directory with raw FASTQs

REFERENCE – reference genome for BLAST

SPECIES – species name

THREADS – number of threads per rule

KMER_LENGTH – k-mer size

MAF – minor allele frequency filter

MAC – minimum allele count

BLOCK_SIZE – for GWAS

KMERS_GWAS_BIN – path to k-mer GWAS binaries

TOOLS_DIR – path to strand annotation binary

PHENO_DIR – directory containing phenotype files

RUN_* – toggle flags for modules (True/False)

Usage
# Activate environment and run Snakemake
snakemake --snakefile master.smk --cores <threads> --configfile config/config.yaml


Ensure all FASTQs are in the correct SOURCE_DIRECTORY.

Run step 07 (extract significant k-mers) manually before the BLAST step.

Final plots and summary tables are generated in plots/.

Notes

The pipeline is modular: each step can be enabled or disabled via RUN_* flags.

Designed for single-species studies at present; multi-species would require minor adjustments.

BLAST step assumes the reference genome has been indexed via makeblastdb.

All intermediate files are organized per sample and per species for reproducibility.

I can also create a diagram/flowchart for the README showing how data moves from raw FASTQs → k-mer counts → GWAS → top k-mers → assembly → BLAST → plots, if you want.

Do you want me to do that next?

do that and one small correction, the part on terminal is done after kmer-gwas

Perfect — thanks for the clarification. I’ve updated the README to reflect that the “extract significant k-mers” terminal step comes after the k-mer GWAS step, not before BLAST. I’ve also added a nice ASCII-style flowchart showing the data flow. Here’s the polished version:

K-mer Based Sex Determination Pipeline – README
Overview

This pipeline performs a k-mer based genome-wide association study (GWAS) to identify k-mers associated with sex in a given species. The workflow includes k-mer counting, strand annotation, kinship correction, GWAS association testing, extraction of significant k-mers, assembly/BLAST of top k-mers, and final visualization.

It is modular and controlled via a Snakemake master workflow, with toggleable steps for flexibility.

Project Structure
PROJECT_DIR/
├── config/
│   └── config.yaml              # Main configuration file
├── workflow/
│   ├── rules/
│   │   ├── 01_files_prep.smk
│   │   ├── 02_kmc.smk
│   │   ├── 03_strand_info.smk
│   │   ├── 04_list_kmers.smk
│   │   ├── 05_kinship_table.smk
│   │   ├── 06_kmer_gwas.smk
│   │   ├── 07_blast_kmers.smk
│   │   └── 09_plot_results.smk
│   └── scripts/
│       ├── plink_to_abyss_kmers.py
│       └── run_plot_kmer_blast.R
├── {SPECIES}/                    # Per-species working directory
│   ├── {sample}/                 # Individual samples
│   │   ├── *_1.fq.gz
│   │   ├── *_2.fq.gz
│   │   ├── output_kmc_*.kmc_pre
│   │   └── kmers_with_strand
│   ├── kmers_list_paths_final.txt
│   ├── kmers_to_use
│   ├── kmers_table.table
│   ├── kmers_table.kinship
│   ├── kmers.assoc
│   ├── kmers.assoc.tab
│   ├── most_assoc_kmers.list
│   └── blast/
│       ├── kmers_blast.out
│       └── mean_contig_len
└── plots/
    └── {SPECIES}_kmer_counts_per_chr.{png,pdf,svg}

Pipeline Modules
01 – Files Preparation (01_files_prep.smk)

Copies raw FASTQ files from the source directory to the project directory.

Handles paired-end samples automatically.

Output: {sample}_1.fq.gz and {sample}_2.fq.gz in the species directory.

02 – K-mer Counting (02_kmc.smk)

Uses KMC to count k-mers in all samples.

Generates both canonical and all k-mers counts.

Output: output_kmc_canon.kmc_pre, output_kmc_all.kmc_pre.

03 – Add Strand Information (03_strand_info.smk)

Adds strand information to k-mer counts.

Uses kmers_add_strand_information binary.

Output: kmers_with_strand.

04 – List K-mers (04_list_kmers.smk)

Combines k-mers from all samples.

Applies MAC (minimum allele count) filtering and removes low-frequency k-mers.

Outputs:

kmers_list_paths_final.txt – paths to per-sample k-mers

dirlist.txt – list of sample directories

kmers_to_use – list of filtered k-mers for GWAS

05 – Kinship Table (05_kinship_table.smk)

Builds a presence/absence k-mer table for all samples.

Computes kinship matrix to correct for population structure using emma_kinship_kmers.

Outputs:

kmers_table.table

kmers_table.kinship

06 – K-mer GWAS (06_kmer_gwas.smk)

Performs GWAS using PLINK with kinship correction.

Converts k-mer table to PLINK format using kmers_table_to_bed.

Outputs:

kmers.assoc – PLINK association results

kmers.assoc.tab – tab-delimited, human-readable table

07 – Extract Significant K-mers (Terminal Step)

After GWAS, manually extract the most significant k-mers by p-value.

Commands:

# Find smallest p-value
cut -f 9 kmers.assoc.tab | grep -v 'P' | sort -g | head -1

# Filter top k-mers (e.g., p <= 3.414e-07)
awk '$9 <= 3.414e-07' kmers.assoc.tab > most_significant_assoc.tab

# Create a list of top k-mers for downstream steps
cat most_significant_assoc.tab | cut -f 2 > most_assoc_kmers.list


Output:

most_significant_assoc.tab

most_assoc_kmers.list

08 – BLAST Top K-mers (07_blast_kmers.smk)

Assembles k-mers using ABYSS.

Maps assembled sequences to reference genome with BLAST.

Outputs:

assoc_kmers_presab.txt

plink_abyss_input.txt

plink_abyss_output.txt

blast/kmers_blast.out

blast/mean_contig_len

09 – Plot Results (09_plot_results.smk)

Generates chromosome-wide visualizations of top k-mers.

Uses R script run_plot_kmer_blast.R.

Outputs:

{SPECIES}_kmer_counts_per_chr.png

{SPECIES}_kmer_counts_per_chr.pdf

{SPECIES}_kmer_counts_per_chr.svg

plots/{SPECIES}_kmer_blast_mapping_summary.csv

Pipeline Flowchart
Raw FASTQs
   │
   ▼
01_files_prep
   │
   ▼
02_kmc (k-mer counting)
   │
   ▼
03_strand_info (add strand info)
   │
   ▼
04_list_kmers (combine & filter k-mers)
   │
   ▼
05_kinship_table (build table & kinship)
   │
   ▼
06_kmer_gwas (GWAS association)
   │
   ▼
07_extract_significant_kmers (terminal)
   │
   ▼
08_blast_kmers (assemble & map)
   │
   ▼
09_plot_results (visualization & summary)

Configuration (config/config.yaml)

Key parameters include:

PROJECT_DIR – base working directory

SOURCE_DIRECTORY – directory with raw FASTQs

REFERENCE – reference genome for BLAST

SPECIES – species name

THREADS – number of threads per rule

KMER_LENGTH – k-mer size

MAF – minor allele frequency filter

MAC – minimum allele count

BLOCK_SIZE – for GWAS

KMERS_GWAS_BIN – path to k-mer GWAS binaries

TOOLS_DIR – path to strand annotation binary

PHENO_DIR – directory containing phenotype files

RUN_* – toggle flags for modules (True/False)

Usage
# Activate environment and run Snakemake
snakemake --snakefile master.smk --cores <threads> --configfile config/config.yaml


Ensure all FASTQs are in the correct SOURCE_DIRECTORY.

Run step 07 (extract significant k-mers) manually after GWAS.

Final plots and summary tables are generated in plots/.

Notes

The pipeline is modular: each step can be enabled or disabled via RUN_* flags.

Designed for single-species studies at present; multi-species would require minor adjustments.

BLAST step assumes the reference genome has been indexed via makeblastdb.

All intermediate files are organized per sample and per species for reproducibility.