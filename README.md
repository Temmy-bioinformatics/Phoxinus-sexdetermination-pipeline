# Sex Determination Analysis Workflows

This repository contains multiple complementary analysis pipelines and scripts focused on sex determination and genomic analyses in *Phoxinus* species.

## Folder Overview

- **coverage-based-analyses/**  
  Workflows and scripts analyzing sequencing coverage differences to identify sex-linked regions.

- **downstream_analysis/**  
  Post-processing, visualization, and integrative analyses of outputs from other pipelines.

- **k-mer-based-analyses/**  
  Pipelines and tools leveraging k-mer based methods for detecting sex-specific sequences and genomic signatures.

- **mitogenome_analyses/**  
  Assembly, annotation, and phylogenetic analyses of mitochondrial genomes.

- **qc-mapping/**  
  Quality control and read mapping workflows, including trimming, alignment, and mapping statistics.

- **snp_calling/**  
  Variant calling pipelines from raw sequencing data, with filtering and quality assessments.

- **snp-based-analyses/**  
  Population genetics, association studies, and other downstream analyses based on SNP data.

## Usage

Each folder contains its own `README.md` with detailed instructions on installation, dependencies, inputs, and running the workflows. Please refer to these for specific guidance.

## Getting Started

- Install necessary tools: Snakemake, samtools, bcftools, R and relevant packages, etc.
- Adjust configuration files in each module as needed.
- Execute workflows using Snakemake or provided scripts.

## Citation

All scripts adapted to snakemake from the original SexFindR pipeline:
Grayson, Phil, et al. "SexFindR: A computational workflow to identify young and old sex chromosomes." BioRxiv (2022): 2022-02.
https://www.biorxiv.org/content/10.1101/2022.02.21.481346v1.full.

