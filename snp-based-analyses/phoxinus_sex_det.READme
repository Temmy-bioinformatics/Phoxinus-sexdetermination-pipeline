---

# Sex Determination Genomic Analyses Pipeline

This repository contains reproducible Snakemake workflows for investigating sex-linked genomic differentiation, nucleotide diversity, SNP density, k-mer-based GWAS, and associated plotting in *Phoxinus* species. The workflows leverage filtered BCF files, raw read data, sex-separated sample lists, and automated Conda environments to ensure reproducibility and scalability.

---

## Directory Structure

```
SexDetermination_snakemake/
├── config/
│   ├── fst_config.yaml
│   ├── sex_nuc_div_config.yaml
│   ├── sex_snp_den_config.yaml
│   ├── snp_gwas_config.yaml
│   ├── kmer_gwas_config.yaml         # Config for all k-mer GWAS steps
│   └── plot_figures_config.yaml      # Config for plotting module
├── resources/
│   ├── sample_lists/
│   │   ├── <species>_Male.vcflist
│   │   └── <species>_Female.vcflist
│   └── bcf_files/
│       └── <species>/Filtered/<species>_allChromosomes.filtered.norm.bcf
├── results/
│   ├── fst/
│   ├── nucleotide_diversity/
│   ├── individual_bcf/
│   ├── snpdensity/
│   ├── SNP_Gwas/
│   │   ├── Bed_files/
│   │   ├── GWAS/
│   │   └── pheno_files/
│   ├── kmer_gwas/
│   │   └── <species>/
│   ├── gwas/
│   └── plotting/                     # Optional folder for storing final plots
├── workflow/
│   ├── envs/
│   │   ├── fst_env.yaml
│   │   ├── p_env.yaml
│   │   ├── bcf_snpden.yaml
│   │   ├── gemma_env.yaml
│   │   ├── kmer_gwas_envs.yaml
│   │   └── r_environment.yaml        # Conda env for R plotting scripts
│   ├── rules/
│   │   ├── fst.smk
│   │   ├── pi.smk
│   │   ├── snp_density.smk
│   │   ├── snp_gwas.smk
│   │   ├── kmer-prep.smk
│   │   ├── kmer-gwas.smk
│   │   ├── top-kmers-analysis.smk
│   │   └── plot_figures.smk          # Plotting module
│   └── Snakefile
└── scripts/
    ├── plot_fst.r
    ├── plot_gwas.r
    ├── plot_snpdensity.r
    ├── plot_nucleotide_diversity.r
    ├── count_windows.py
    └── others...
```

---

## Configuration Files

* `kmer_gwas_config.yaml`: Settings for all k-mer GWAS steps
* `plot_figures_config.yaml`: Paths, scripts, and parameters for all plotting workflows
* Other config files for FST, nucleotide diversity, SNP density, and SNP GWAS as before

---

## Conda Environments

* Separate environments under `workflow/envs/` for modular reproducibility:

  * FST, nucleotide diversity, SNP density, SNP GWAS, k-mer GWAS, and plotting (R)
* The plotting environment (`r_environment.yaml`) includes all R packages for figure generation (ggplot2, patchwork, etc.)
* Environments are automatically activated by Snakemake’s `--use-conda` flag

---

## Pipeline Modules Overview

### K-mer GWAS Pipelines

* **kmer-prep.smk:** FASTQ copying, k-mer counting, strand info
* **kmer-gwas.smk:** K-mer table building, kinship, GWAS via GEMMA
* **top-kmers-analysis.smk:** Filtering top hits, contig assembly, BLAST annotation

### Population Genomic Metrics

* FST, nucleotide diversity, SNP density, SNP-based GWAS pipelines as before

### Plotting Module (`plot_figures.smk`)

* Generates high-quality publication-ready plots and summary text files for multiple analyses:

  * **FST**: Windowed FST values plotted along chromosomes with top 5% and 1% cutoffs, including thresholds saved as text files.
  * **GWAS**: Manhattan-style plots of SNP association p-values, with top SNP subsets exported for downstream analysis.
  * **SNP Density**: Plots of variant density differences by sex, plus proportion plots of significant SNP density windows.
  * **Nucleotide Diversity**: Scatterplots showing male-female π differences and variant count differences per 10kb window.

* Each plotting script outputs:

  * **TIFF/PDF figures** for visualization
  * **Text summary files** (e.g., thresholds, significant SNP lists, merged data tables) for reproducibility and downstream integration

---

## Output Files and SexFindR Integration

The plotting workflows produce **standardized text output files** that serve as **input for the general SexFindR plotting scripts**. These include, but are not limited to:

| File Example                                                   | Description                                   | Use Case                                         |
| -------------------------------------------------------------- | --------------------------------------------- | ------------------------------------------------ |
| `<species>/<species>_Fst_10k_thresholds.txt`                   | Numeric FST cutoff thresholds (5% and 1%)     | For highlighting significant genomic regions     |
| `<species>/<species>_NucDiv_10kb_sex_differences.txt`          | Nucleotide diversity differences by sex       | For combined nucleotide diversity plots          |
| `<species>/plot/<species>_SNPdensity_SexFindR.txt`             | SNP density windows with sex-separated values | Input for SNP density visualization              |
| `<species>/plot/<species>_significant_SNPdensity_SexFindR.txt` | Significant SNP density windows (p ≤ 0.001)   | Highlighting strongly sex-linked genomic windows |
| `<species>/output/<species>_gwas_plot_values.txt`              | GWAS SNP association p-values and positions   | Used for integrated GWAS signal plotting         |

These files allow integrated visualization of overlapping signals from diverse population genomic metrics, facilitating robust sex-linked region discovery.

---

## Running the Full Workflow

From the project root, run the entire pipeline:

```bash
snakemake --cores 10 --use-conda
```

To run individual modules selectively (e.g., k-mer GWAS or plotting):

```bash
snakemake -s workflow/rules/kmer-prep.smk --cores 10 --use-conda
snakemake -s workflow/rules/kmer-gwas.smk --cores 10 --use-conda
snakemake -s workflow/rules/top-kmers-analysis.smk --cores 10 --use-conda
snakemake -s workflow/rules/plot_figures.smk --cores 4 --use-conda
```

---

## Recommendations and Next Steps

* Use the standardized text files generated by plotting workflows as inputs for the integrated SexFindR plotting scripts that combine multiple signals.
* Visualize candidate sex-linked regions with multi-track genome browser views or composite plots.
* Extend analyses to additional populations or species for comparative evolutionary insights.

---
