Here's a detailed README file for your Snakemake master workflow based on the `master.smk` you shared. This README explains the purpose, setup, configuration, usage, and output expectations to guide users through the workflow.

---

# Sex Determination k-mer GWAS Snakemake Workflow

## Overview

This Snakemake workflow performs a k-mer-based Genome-Wide Association Study (GWAS) focused on sex determination in the *Phoxinus* genus (or other species of interest). The pipeline is modular and controlled via a single configuration file, allowing users to toggle individual analysis steps such as k-mer counting, kinship matrix calculation, association testing, and BLAST annotation.

The workflow is designed to handle Illumina sequencing FASTQ files organized by species and gender, and supports multi-threaded execution to efficiently process large datasets.

---

## Features

* Modular pipeline with independent toggles for each analysis step
* Supports multiple species and sample genders
* Uses a reference genome for kmers mapping and annotation
* Performs k-mer counting with KMC or similar tools
* Calculates kinship matrix to account for relatedness among samples
* Performs k-mer based GWAS for association with sex or other phenotypes
* BLAST of significant kmers for genomic localization and annotation
* Automatic sample discovery from structured directories

---

## Repository Structure (example)

```text
.
├── config/
│   └── config.yaml              # Main configuration file controlling paths and parameters
├── workflow/
│   └── rules/
│       ├── 01_files_prep.smk
│       ├── 02_kmc.smk
│       ├── 03_strand_info.smk
│       ├── 04_list_kmers.smk
│       ├── 05_kinship_table.smk
│       ├── 06_kmer_gwas.smk
│       └── 07_blast_kmers.smk
├── master.smk                   # Master Snakefile that orchestrates workflow based on config toggles
└── README.md                    # This documentation file
```

---

## Prerequisites

* **Snakemake** (tested with Snakemake ≥5.0)
* Python 3 environment
* Required bioinformatics tools installed and accessible (e.g., KMC, BLAST, kmerGWAS binaries)
* Reference genome FASTA file (rearranged as needed)
* Phenotype files in the appropriate format
* Sequencing data organized by species and gender

---

## Configuration (`config.yaml`)

The workflow behavior is controlled by the `config.yaml` file, located (by default) at:

```
/home/toriowo/Snakemake_projects/sexfindr-snakemake/k-mer-based-analyses/config/config.yaml
```

### Important configuration fields:

| Parameter          | Description                                                                     |
| ------------------ | ------------------------------------------------------------------------------- |
| `SOURCE_DIRECTORY` | Path to directory containing raw FASTQ files organized by species/gender/sample |
| `PROJECT_DIR`      | Base directory to write all pipeline output files                               |
| `TOOLS_DIR`        | Directory containing kmerGWAS binaries and scripts                              |
| `PHENO_DIR`        | Directory containing phenotype files                                            |
| `REFERENCE`        | Path to reference genome FASTA                                                  |
| `KMERS_GWAS_BIN`   | Path to kmerGWAS executables                                                    |
| `SPECIES`          | Single species name or list of species to analyze                               |
| `THREADS`          | Number of CPU threads to use                                                    |
| `KMER_LENGTH`      | Length of kmers for analysis (e.g., 31)                                         |
| `MAC`              | Minimum allele count for kmers to retain                                        |
| `MAF`              | Minor allele frequency filter                                                   |
| `BLOCK_SIZE`       | Block size for splitting genome or data                                         |

### Analysis toggles (set `true` or `false`):

* `RUN_FILES_PREP`: Prepare input files (e.g., filtering)
* `RUN_KMC`: Run kmer counting
* `RUN_STRAND_INFO`: Extract strand bias information for kmers
* `RUN_LIST_KMERS`: Generate kmer lists for downstream analysis
* `RUN_KINSHIP_TABLE`: Calculate kinship/relatedness matrix
* `RUN_KMER_GWAS`: Perform the kmer GWAS analysis
* `RUN_BLAST_KMERS`: Run BLAST to annotate significant kmers

---

## Directory Structure of Input Data

The workflow expects raw FASTQ files to be organized as:

```
<SOURCE_DIRECTORY>/<SPECIES>/<Gender>/<sample>_1.fq.gz
<SOURCE_DIRECTORY>/<SPECIES>/<Gender>/<sample>_2.fq.gz
```

Example:

```
/share/pool/toriowo/Resources/Sex_Determination/Combined_Metadata/Fastq_Files_List/RhinePhoxinus/Male/Sample01_1.fq.gz
/share/pool/toriowo/Resources/Sex_Determination/Combined_Metadata/Fastq_Files_List/RhinePhoxinus/Male/Sample01_2.fq.gz
/share/pool/toriowo/Resources/Sex_Determination/Combined_Metadata/Fastq_Files_List/RhinePhoxinus/Female/Sample02_1.fq.gz
```

The Snakefile will dynamically detect samples based on this directory structure.

---

## How to Run the Workflow

1. **Edit `config.yaml`**
   Customize paths, parameters, and toggle the analysis steps you want to run.

2. **Launch Snakemake**

```bash
snakemake --snakefile master.smk --cores <number_of_threads> --use-conda
```

Replace `<number_of_threads>` with the number of CPU cores you want to allocate.

3. **Workflow Behavior**

* The master Snakefile imports configuration and dynamically includes Snakefiles for each step based on toggles.
* Final outputs per enabled module are aggregated to the `all` rule.
* Progress and executed modules are printed to `stderr` for easy monitoring.

---

## Output

Outputs are organized under:

```
<PROJECT_DIR>/<SPECIES>/<sample>/
```

Examples of outputs generated depending on enabled steps:

* Preprocessed FASTQ files
* KMC kmer count databases (`.kmc_pre`)
* Kmer strand information files
* Lists of kmers for GWAS
* Kinship matrices (`.table` and `.kinship`)
* GWAS association result files (`kmers.assoc`, `kmers.assoc.tab`)
* BLAST annotation results for top kmers

---

## Extending or Modifying the Workflow

* To add or modify steps, create or edit rule files in the `workflow/rules/` directory.
* Toggle your desired steps in `config.yaml`.
* Adjust sample discovery or reference paths in `master.smk` if your data layout differs.
* Consider adding logging or additional QC steps as needed.

---

