# SNP-Based Sex Determination Analysis Pipeline

A modular Snakemake pipeline for identifying sex-linked genomic regions through comprehensive population genomic analyses including FST, nucleotide diversity, SNP density, and GWAS.

> **Status:** Production-ready | Fully tested | Actively maintained

---

## 📋 Table of Contents

- [Overview](#overview)
- [Quick Start](#quick-start)
- [Pipeline Modules](#pipeline-modules)
- [Directory Structure](#directory-structure)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [Output Files](#output-files)
- [Workflow Details](#workflow-details)
- [Species Customization](#species-customization)
- [Troubleshooting](#troubleshooting)
- [Performance](#performance)
- [Citation](#citation)

---

## 🎯 Overview

This pipeline identifies candidate sex-linked regions by comparing male and female genomic patterns across multiple statistical approaches. It is designed for non-model organisms where sex chromosomes are unknown or poorly characterized.

### Key Features

- ✅ **Modular Design** - Run individual analyses or the complete pipeline
- ✅ **Automated Environment Management** - Conda handles all dependencies
- ✅ **Species-Agnostic** - Customizable color schemes and scalable to multiple species
- ✅ **Publication-Ready** - Generates high-quality plots and publication tables
- ✅ **Integrated Analysis** - Identifies regions significant across multiple methods
- ✅ **Reproducible** - Version-controlled environments and workflows

### Analyses Performed

| Module | Analysis | Output |
|--------|----------|--------|
| **01-FST** | Fixation index between sexes | Manhattan plots, candidate SNPs, window enrichment |
| **02-Pi** | Nucleotide diversity differences | Sex-difference plots, diversity metrics |
| **03-BCF** | Individual sample extraction | Per-sample BCF files for downstream analysis |
| **04-SNP Density** | Differential variant density | Permutation-tested significant windows |
| **05-GWAS** | Sex association analysis | Association statistics, Manhattan plots |
| **06-Integration** | Cross-method validation | Candidate regions, correlations, combined plots |

---

## 🚀 Quick Start

```bash
# 1. Clone repository
cd /path/to/your/project

# 2. Install Snakemake
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake

# 3. Configure pipeline
cp config/config.yaml.example config/config.yaml
# Edit config.yaml with your paths

# 4. Dry run to test
snakemake --use-conda --cores 8 -n

# 5. Run complete pipeline
snakemake --use-conda --cores 8

# 6. Or run specific modules
snakemake --use-conda --cores 8 --config RUN_FST=True RUN_NUC_DIV=True
```

---

## 📦 Pipeline Modules

### Module 1: FST Analysis (`01-fst.smk`)

**Purpose:** Identifies regions with high genetic differentiation between males and females.

**Workflow:**
1. Calculate per-SNP FST (Weir & Cockerham)
2. Calculate FST in 10kb sliding windows
3. Extract top 5% FST SNPs
4. Sort SNPs by genomic coordinates
5. Count high-FST SNPs per window
6. Generate Manhattan plots with 95th/99th percentile thresholds

**Key Outputs:**
- `{species}_Fst.txt.weir.fst` - Per-SNP FST values
- `{species}_Fst_10k.windowed.weir.fst` - Windowed FST
- `{species}_Fst_10k_plot.tiff` - Manhattan plot
- `{species}_Fst_window_count.txt` - Window enrichment counts

**Statistical Approach:**
- Weir & Cockerham's FST estimator
- Thresholds at 95th and 99th percentiles
- Window-based enrichment for identifying candidate regions

---

### Module 2: Nucleotide Diversity (`02-pi.smk`)

**Purpose:** Compares genetic diversity between sexes to identify regions with sex-specific patterns (e.g., reduced diversity on sex chromosomes in the heterogametic sex).

**Workflow:**
1. Calculate π in 10kb windows for males
2. Calculate π in 10kb windows for females
3. Compute male - female differences
4. Generate genome-wide difference plots

**Key Outputs:**
- `{species}_Male_10k_nucleotide_diversity.windowed.pi` - Male diversity
- `{species}_Female_10k_nucleotide_diversity.windowed.pi` - Female diversity
- `{species}_NucDiv_10kb_sex_differences.txt` - Difference table
- `{species}_pi_10kb_window.tiff` - Nucleotide diversity plot
- `{species}_var_per_10kb_window.tiff` - Variant count plot

**Statistical Approach:**
- Watterson's θ estimator
- Direct comparison of male vs female diversity
- Identifies regions with reduced diversity in one sex

---

### Module 3: BCF Extraction (`03-extract-sample-bcf.smk`)

**Purpose:** Extracts individual sample genotypes from multi-sample BCF for per-sample analyses.

**Workflow:**
1. Extract individual BCF per sample
2. Index BCF files with bcftools

**Key Outputs:**
- `{species}/{gender}/{sample}.bcf` - Individual sample BCF
- `{species}/{gender}/{sample}.bcf.csi` - BCF index

**Note:** This is a prerequisite for SNP density analysis.

---

### Module 4: SNP Density Analysis (`04-run-snp-density.smk`)

**Purpose:** Identifies genomic windows with significantly different SNP densities between sexes using rigorous permutation testing.

**Workflow:**
1. Calculate SNP density in 10kb windows per sample
2. Compute mean density for males and females
3. Perform 1000 permutation tests
4. Calculate empirical p-values
5. Identify significant windows (p ≤ 0.001)
6. Plot proportion of significant SNPs per chromosome

**Key Outputs:**
- `{gender}_{sample}_10k.snpden` - Per-sample SNP density
- `{species}_SNPdensity_SexFindR.txt` - All windows with p-values
- `{species}_significant_SNPdensity_SexFindR.txt` - Significant windows only
- `{species}_proportion_p001.pdf` - Significance plot

**Statistical Approach:**
- 1000 permutation iterations
- Empirical p-value calculation
- FDR correction at p ≤ 0.001
- Controls for genome-wide variation patterns

---

### Module 5: SNP GWAS (`05-snp-gwas.smk`)

**Purpose:** Performs genome-wide association analysis treating sex as a binary phenotype.

**Workflow:**
1. Create phenotype file (males=1, females=2)
2. Convert BCF to PLINK format
3. Run GEMMA linear mixed model
4. Extract top 5% associated SNPs
5. Sort and count SNPs in windows
6. Generate Manhattan plots

**Key Outputs:**
- `{species}.pheno.txt` - Phenotype file
- `{species}.bed/.bim/.fam` - PLINK binary format
- `{species}_gemma.assoc.txt` - Association statistics
- `{species}_GWAS.tiff` - Manhattan plot
- `{species}_GWAS_window_count.txt` - Window enrichment

**Statistical Approach:**
- GEMMA linear mixed model
- Controls for population structure
- Identifies SNPs associated with sex determination
- Window-based enrichment analysis

---

### Module 6: Integrated Analysis (`06-plot-results.smk`)

**Purpose:** Integrates results from all four analyses to identify high-confidence candidate regions and assess method concordance.

**Workflow:**
1. Rank windows from each analysis method
2. Calculate Spearman correlations between methods
3. Identify windows in top 100 across multiple methods
4. Generate 4-panel integrated plot
5. Export candidate region tables

**Key Outputs:**
- `candidate_windows_SNP_GWAS_Fst_pi.txt` - Top 100 in all 4 methods
- `candidate_windows_GWAS_Pi_Fst.txt` - Top 100 in 3 methods (GWAS+Pi+FST)
- `candidate_windows_SNP_Fst_Pi.txt` - Top 100 in 3 methods (SNP+FST+Pi)
- `{species}_correlation_rankings.txt` - Method correlations
- `{species}_integrated_plot.png` - 4-panel combined visualization

**Statistical Approach:**
- Ranking-based integration
- Spearman correlation for method concordance
- Multiple filtering strategies for candidate selection
- Visual integration with shared genomic coordinates

---

## 📂 Directory Structure

```
snp-based-analyses/
├── Snakefile                              # Master workflow controller
├── config/
│   └── config.yaml                        # Main configuration file
├── workflow/
│   ├── rules/
│   │   ├── 01-fst.smk                    # FST analysis
│   │   ├── 02-pi.smk                     # Nucleotide diversity
│   │   ├── 03-extract-sample-bcf.smk     # BCF extraction
│   │   ├── 04-run-snp-density.smk        # SNP density
│   │   ├── 05-snp-gwas.smk               # GWAS analysis
│   │   └── 06-plot-results.smk           # Integrated plotting
│   ├── envs/
│   │   ├── fst_env.yaml                  # VCFtools environment
│   │   ├── r-plots_env.yaml              # R plotting environment
│   │   ├── gemma_env.yaml                # GEMMA/PLINK environment
│   │   └── p_env.yaml                    # Python environment
│   └── scripts/
│       ├── r-scripts/
│       │   ├── plot_fst.r                # FST plotting
│       │   ├── extract_top_fst.R         # Extract top FST SNPs
│       │   ├── plot_nucleotide_diversity.R  # Pi plotting
│       │   ├── plot_snp_density.R        # SNP density plotting
│       │   ├── plot_snp-gwas.R           # GWAS plotting
│       │   └── sexfindr-combined-plot.R  # Integrated plotting
│       └── count_windows.py              # Window counting utility
├── resources/
│   ├── sample_lists/
│   │   ├── {species}_Male.vcflist        # Male sample IDs
│   │   └── {species}_Female.vcflist      # Female sample IDs
│   ├── bcf_files/
│   │   └── {species}/Filtered/{species}_allChromosomes.filtered.norm.bcf
│   └── reference/
│       └── {reference}.fasta.fai         # Reference genome index
└── results/
    ├── fst/{species}/                    # FST results
    ├── pi/{species}/                     # Nucleotide diversity results
    ├── individual_bcf/{species}/{gender}/  # Individual BCF files
    ├── snpdensity/{species}/             # SNP density results
    ├── snp_gwas/                         # GWAS results
    │   ├── pheno_files/{species}/
    │   ├── bed_files/{species}/
    │   └── gwas/{species}/output/
    ├── plotting/{species}/               # Individual plots
    └── combined_plots/{species}/         # Integrated analysis
```

---

## 🔧 Installation

### Prerequisites

- **Operating System:** Linux or macOS
- **Conda/Mamba:** Miniforge recommended (https://github.com/conda-forge/miniforge)
- **Snakemake:** ≥ 7.0
- **Disk Space:** ~5-10 GB per species (varies with genome size)
- **CPU:** Multi-core system recommended (pipeline is highly parallelizable)
- **Memory:** 8-16 GB RAM minimum

### Setup Steps

```bash
# 1. Install Miniforge (if not already installed)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh

# 2. Create Snakemake environment
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake

# 3. Verify installation
snakemake --version

# 4. Prepare your data
# - Multi-sample BCF file (filtered and normalized)
# - Sample lists (one sample ID per line)
# - Reference genome index (.fai file)
```

---

## ⚙️ Configuration

### Main Configuration File: `config/config.yaml`

```yaml
#####################################
# PROJECT PATHS
#####################################
PROJECT_DIR: "/path/to/project"
OUTPUT_DIR: "results"
BCF_DIR: "resources/bcf_files"
SAMPLE_DIR: "resources/sample_lists"
INDV_BCF_DIR: "results/individual_bcf"

#####################################
# SPECIES
#####################################
SPECIES:
  - SeinePhoxinus
  - RhinePhoxinus
  # Add more species as needed

#####################################
# SAMPLE LISTS
#####################################
VCF_LISTS:
  - SeinePhoxinus_Male.vcflist
  - SeinePhoxinus_Female.vcflist
  - RhinePhoxinus_Male.vcflist
  - RhinePhoxinus_Female.vcflist

#####################################
# REFERENCE FILES
#####################################
ref_file: "resources/reference/Final_Phoxy_1_rearranged.fasta.fai"

#####################################
# R SCRIPTS
#####################################
fst_rscript: "workflow/scripts/r-scripts/plot_fst.r"
extract_fst_rscript: "workflow/scripts/r-scripts/extract_top_fst.R"
nuc_diversity_rscript: "workflow/scripts/r-scripts/plot_nucleotide_diversity.R"
snp_density_rscript: "workflow/scripts/r-scripts/plot_snp_density.R"
gwas_rscript: "workflow/scripts/r-scripts/plot_snp-gwas.R"
integrated_analysis_rscript: "workflow/scripts/r-scripts/sexfindr-combined-plot.R"

#####################################
# PYTHON UTILITIES
#####################################
count_window: "scripts/count_windows.py"

#####################################
# RESOURCES
#####################################
THREADS: 8

#####################################
# MODULE FLAGS
#####################################
RUN_FST: True
RUN_NUC_DIV: True
RUN_EXTRACT_BCF: True
RUN_SNP_DENSITY: True
RUN_SNP_GWAS: True
RUN_PLOTTING: True
```

### Input File Formats

**Sample Lists (`.vcflist`):**
```
Sample001
Sample002
Sample003
```
One sample ID per line, matching exactly the sample names in your BCF file.

**BCF Files:**
- Multi-sample BCF format
- Filtered and normalized
- Must include all samples (males and females)
- Must have corresponding index (.csi)

**Reference Index (.fai):**
```
Chr1    45812323    52    60    61
Chr2    32451234    45912435    60    61
Chr3    28934567    78367669    60    61
```
Standard samtools faidx format.

---

## 🎮 Usage

### Complete Pipeline

```bash
# Activate environment
conda activate snakemake

# Run all enabled modules
snakemake --use-conda --cores 8

# Dry run to preview execution
snakemake --use-conda --cores 8 -n

# Generate workflow diagram
snakemake --dag | dot -Tpdf > workflow.pdf
```

### Module-Specific Execution

```bash
# FST analysis only
snakemake --use-conda --cores 8 \
  --config RUN_FST=True RUN_NUC_DIV=False RUN_EXTRACT_BCF=False \
  RUN_SNP_DENSITY=False RUN_SNP_GWAS=False RUN_PLOTTING=False

# FST and Nucleotide Diversity
snakemake --use-conda --cores 8 \
  --config RUN_FST=True RUN_NUC_DIV=True

# SNP Density and GWAS (requires BCF extraction)
snakemake --use-conda --cores 8 \
  --config RUN_EXTRACT_BCF=True RUN_SNP_DENSITY=True RUN_SNP_GWAS=True

# Integrated analysis (requires all modules)
snakemake --use-conda --cores 8 --config RUN_PLOTTING=True
```

### HPC/Cluster Execution

For SLURM-based clusters:

```bash
snakemake --use-conda --cores 100 \
  --cluster "sbatch -p {cluster.partition} -t {cluster.time} \
  -c {threads} --mem={resources.mem_mb}" \
  --cluster-config cluster.yaml \
  --jobs 50
```

Example `cluster.yaml`:
```yaml
__default__:
  partition: medium
  time: "24:00:00"
  mem_mb: 8000

run_gemma:
  partition: high
  time: "48:00:00"
  mem_mb: 16000
```

---

## 📊 Output Files

### FST Module (`results/fst/{species}/`)

```
{species}_Fst.txt.weir.fst                    # Per-SNP FST values
{species}_Fst_10k.windowed.weir.fst           # 10kb windowed FST
{species}__all_Fst_sort_cut.txt               # Top 5% FST SNPs (unsorted)
{species}_coordinated_sorted_Fst_cut.txt      # Top 5% FST SNPs (sorted)
{species}_Fst_window_count.txt                # Window enrichment counts
```

```
results/plotting/{species}/
{species}_Fst_10k_plot.tiff                   # Manhattan plot
{species}_Fst_10k_thresholds.txt              # 95th/99th percentile values
```

### Nucleotide Diversity Module (`results/pi/{species}/`)

```
{species}_Male_10k_nucleotide_diversity.windowed.pi    # Male π values
{species}_Female_10k_nucleotide_diversity.windowed.pi  # Female π values
{species}_NucDiv_10kb_sex_differences.txt              # Difference table
{species}_pi_10kb_window.tiff                          # π difference plot
{species}_var_per_10kb_window.tiff                     # Variant count plot
```

### SNP Density Module (`results/snpdensity/{species}/`)

```
Male_{sample}_10k.snpden                      # Per-sample density (males)
Female_{sample}_10k.snpden                    # Per-sample density (females)

plot/
├── {species}_proportion_p001.pdf             # Significance plot
├── {species}_SNPdensity_SexFindR.txt        # All windows + p-values
└── {species}_significant_SNPdensity_SexFindR.txt  # Significant only
```

### GWAS Module (`results/snp_gwas/`)

```
pheno_files/{species}/{species}.pheno.txt     # Phenotype file

bed_files/{species}/
├── {species}.bed                             # PLINK binary
├── {species}.bim                             # Variant info
└── {species}.fam                             # Sample info

gwas/{species}/output/
├── {species}_gemma.assoc.txt                 # Association results
├── {species}_gemma.log                       # GEMMA log
├── {species}_GWAS.tiff                       # Manhattan plot
├── {species}_gwas_plot_values.txt            # Plot data
├── {species}_GWAS_out_sort_cut.txt           # Top 5% SNPs
├── {species}_coordinated_sorted_GWAS_cut.txt # Sorted top SNPs
└── {species}_GWAS_window_count.txt           # Window enrichment
```

### Integrated Analysis (`results/combined_plots/{species}/`)

```
candidate_windows_SNP_GWAS_Fst_pi.txt         # Top 100 in all 4 methods
candidate_windows_GWAS_Pi_Fst.txt             # Top 100 in GWAS+Pi+FST
candidate_windows_SNP_Fst_Pi.txt              # Top 100 in SNP+FST+Pi
{species}_correlation_rankings.txt            # Spearman correlations
{species}_integrated_plot.png                 # 4-panel combined plot
```

---

## 🔬 Workflow Details

### Module Dependencies

```
BCF File + Sample Lists
    │
    ├─► FST Analysis (independent) ─────────────┐
    │                                            │
    ├─► Nucleotide Diversity (independent) ─────┤
    │                                            │
    ├─► Extract Individual BCF                  │
    │       │                                    │
    │       └─► SNP Density ────────────────────┤
    │                                            │
    └─► SNP GWAS (independent) ─────────────────┤
                                                 │
                                                 ▼
                                        Integrated Analysis
                                    (requires all 4 analyses)
```

### Statistical Approaches

#### FST (Fixation Index)
- **Method:** Weir & Cockerham's FST estimator
- **Interpretation:** Measures allele frequency differences between sexes
- **Thresholds:** 95th percentile (top 5%) and 99th percentile (top 1%)
- **Use Case:** Identifies regions under divergent selection or with restricted gene flow

#### Nucleotide Diversity (π)
- **Method:** Watterson's θ estimator
- **Interpretation:** Measures within-population diversity
- **Analysis:** Male - Female comparison
- **Use Case:** Sex chromosomes often show reduced diversity in heterogametic sex (e.g., Y chromosome in XY systems)

#### SNP Density
- **Method:** Permutation testing (1000 iterations)
- **Interpretation:** Identifies windows with differential variant accumulation
- **Significance:** Empirical p-values ≤ 0.001
- **Use Case:** Regions with reduced recombination often show variant clustering

#### GWAS
- **Method:** GEMMA linear mixed model
- **Interpretation:** Association between genotype and sex phenotype
- **Controls:** Population structure via kinship matrix
- **Use Case:** Direct identification of sex-determining loci

#### Integration
- **Method:** Ranking-based across all methods
- **Correlation:** Spearman's ρ for method concordance
- **Filtering:** Multiple strategies (all 4 methods, 3-way combinations)
- **Use Case:** High-confidence candidate regions

---

## 🎨 Species Customization

### Adding New Species

1. **Update `config.yaml`:**
```yaml
SPECIES:
  - ExistingSpecies
  - NewSpecies  # Add here

VCF_LISTS:
  - ExistingSpecies_Male.vcflist
  - ExistingSpecies_Female.vcflist
  - NewSpecies_Male.vcflist        # Add here
  - NewSpecies_Female.vcflist      # Add here
```

2. **Create sample lists:**
```bash
# Create resources/sample_lists/NewSpecies_Male.vcflist
# Create resources/sample_lists/NewSpecies_Female.vcflist
```

3. **Add species colors to R scripts:**

In each plotting script (`plot_fst.r`, `plot_nucleotide_diversity.R`, `plot_snp_density.R`, `plot_snp-gwas.R`, `sexfindr-combined-plot.R`):

```r
color_map <- list(
    SeinePhoxinus = c("#9897A9", "#ab9519"),
    RhinePhoxinus = c("#9897A9", "#66A61E"),
    NewSpecies = c("#color1", "#color2")  # Add your colors here
)
```

**Color Tips:**
- Use hex colors for consistency
- First color for even chromosomes, second for odd
- Choose contrasting colors for visibility
- Test colors for colorblind accessibility

---

## 🐛 Troubleshooting

### Common Issues

#### 1. Missing Conda Environments
```bash
Error: Environment file not found

Solution:
# Verify environment files exist
ls workflow/envs/

# Should contain:
# - fst_env.yaml
# - r-plots_env.yaml
# - gemma_env.yaml
# - p_env.yaml
```

#### 2. Sample Name Mismatches
```bash
Error: Sample names in VCF don't match sample lists

Solution:
# Check sample names in BCF
bcftools query -l resources/bcf_files/{species}/Filtered/{species}_allChromosomes.filtered.norm.bcf

# Compare with sample lists
cat resources/sample_lists/{species}_Male.vcflist
cat resources/sample_lists/{species}_Female.vcflist

# Names must match exactly (case-sensitive)
```

#### 3. Conda Path Issues
```bash
Error: Conda environment path not found

Solution:
# Conda environment paths in rule files are relative to the rule file location
# From workflow/rules/, paths should be: envs/fst_env.yaml
# NOT: ../envs/fst_env.yaml
# Snakemake resolves paths from the rule file directory
```

#### 4. Memory Issues
```bash
Error: GEMMA runs out of memory

Solution 1: Increase memory allocation
# Edit 05-snp-gwas.smk
resources:
    mem_mb = 16000  # Increase from 8000

Solution 2: Filter variants more stringently
# Use stricter MAF or missingness filters before running pipeline
```

#### 5. Empty Output Files
```bash
Error: Output files created but empty

Solution:
# Check input file quality
bcftools stats input.bcf

# Verify sufficient variants
# Ensure balanced male/female samples
# Check for filtering issues
```

#### 6. R Package Errors
```bash
Error: R package not found

Solution:
# Activate the conda environment
conda activate .snakemake/conda/{env_hash}

# Test package availability
R -e "library(tidyverse)"

# If missing, environment may need rebuilding
rm -rf .snakemake/conda/{env_hash}
snakemake --use-conda --cores 8
```

### Debugging Tips

**Check Snakemake Logs:**
```bash
cat .snakemake/log/{timestamp}.snakemake.log
```

**Check Rule-Specific Logs:**
```bash
# SNP density logs
cat results/snpdensity/{species}/logs/*.log

# GEMMA logs
cat results/snp_gwas/gwas/{species}/output/{species}_gemma.log
```

**Validate Input Files:**
```bash
# Check BCF integrity
bcftools view input.bcf | head

# Count variants
bcftools view -H input.bcf | wc -l

# Verify sample lists
wc -l resources/sample_lists/*.vcflist

# Check reference index
head resources/reference/*.fai
```

**Rerun Failed Jobs:**
```bash
# Rerun incomplete files
snakemake --use-conda --cores 8 --rerun-incomplete

# Force rerun specific rule
snakemake --use-conda --cores 8 --forcerun plot_fst

# Unlock working directory if needed
snakemake --unlock
```


## ⚡ Performance Optimization

### Parallelization

The pipeline is highly parallelizable. Adjust based on your system:


### Resource Allocation

**Recommended resources per module:**

| Module | Cores | Memory | Time |
|--------|-------|--------|------|
| FST | 1-4 | 4-8 GB | 1-4 hours |
| Nucleotide Diversity | 1-4 | 4-8 GB | 1-4 hours |
| BCF Extraction | 1 per sample | 2-4 GB | 10-30 min/sample |
| SNP Density | 1 per sample | 2-4 GB | 10-30 min/sample |
| GWAS | 4-8 | 8-16 GB | 2-6 hours |
| Integrated Plotting | 1 | 4-8 GB | 10-30 min |

### Disk Space Requirements

Estimate **5-10 GB per species** for all results, depending on:
- Genome size
- Number of samples
- Variant density
- Number of analyses run
