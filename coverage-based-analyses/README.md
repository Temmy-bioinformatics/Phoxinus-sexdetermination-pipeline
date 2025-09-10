# BAM Coverage Comparison Pipeline

This repository contains a modular **Snakemake** workflow for BAM-based coverage comparison analyses.
It is designed to run coverage statistics, merge BAMs, compute differential coverage, and optionally run `DNAcopy` for downstream segmentation.

---

## Table of Contents

* [Features](#features)
* [Pipeline Overview](#pipeline-overview)
* [Installation](#installation)
* [Configuration](#configuration)
* [Running the Pipeline](#running-the-pipeline)
* [Modules](#modules)
* [Outputs](#outputs)
* [Notes](#notes)

---

## Features

* Modular **Snakefiles**, each responsible for a distinct analysis step
* Toggle modules on/off via `config.yaml`
* Automates merging, coverage statistics, and differential coverage calculations
* Optional **DNAcopy** segmentation for detecting copy number differences
* Fully compatible with **Conda environments**

---

## Pipeline Overview

The workflow consists of six modules:

1. **Samtools Stats** – Compute coverage and mapping statistics for all BAM files
2. **MultiQC** – Aggregate QC metrics into an HTML report
3. **Merge BAMs** – Merge BAMs by species and sex, indexing and computing stats
4. **Coverage Stats** – Compute modal and min/max coverage values for male and female BAMs
5. **DifCover** – Perform differential coverage calculation between sexes
6. **DNAcopy Input** – Generate log2-adjusted ratios for use in DNAcopy segmentation

---

## Installation

1. Clone this repository:

```bash
git clone https://github.com/yourusername/bam-coverage-pipeline.git
cd bam-coverage-pipeline
```

2. Install **Snakemake** (via Conda recommended):

```bash
conda install -c bioconda snakemake
```

3. Make sure the environments in `envs/` are installed or will be created automatically.

---

## Configuration

Set paths and parameters in `config.yaml`:

```yaml
# Paths
PROJECT_DIR: "/path/to/project/results"
OUTPUT_DIR: "/path/to/output"
MERGED_BAM_DIR: "/path/to/merged_bams"
SCRIPTS_DIR: "/path/to/scripts"

# Parameters
SPECIES: ["RhinePhoxinus", "DanubeCsikii"]
THREADS: 8

# Module toggles
RUN_SAMTOOLS_STATS: true
RUN_MULTIQC: true
RUN_MERGE_BAMS: true
RUN_COV_STATS: true
RUN_DIFCOVER: true
RUN_DNACOPY: true

# Optional R script execution
RUN_R_DNACOPY: true
```

*You can toggle individual modules by setting their `RUN_*` values to `false`.*

---

## Running the Pipeline

Run with:

```bash
snakemake --cores 16 --use-conda -s master.smk
```

The pipeline will only run the modules enabled in your config file.

---

## Modules

| Module Name    | Snakefile                 | Description                               |
| -------------- | ------------------------- | ----------------------------------------- |
| Samtools Stats | `rules/01-samtools.smk`   | Compute BAM statistics                    |
| MultiQC        | `rules/02-multiqc.smk`    | Aggregate QC results into an HTML report  |
| Merge BAMs     | `rules/03-merge-bams.smk` | Merge BAMs by species and sex             |
| Coverage Stats | `rules/04-cov-stats.smk`  | Compute coverage metrics and ratios       |
| DifCover       | `rules/05-difcover.smk`   | Calculate differential coverage           |
| DNAcopy Input  | `rules/06-dnacopy.smk`    | Generate log2-adjusted ratios for DNAcopy |

---

## Optional DNAcopy R Script

The `dnacopy_input` rule generates `.log2adj.txt` files.
If `RUN_R_DNACOPY` is set to `true` in `config.yaml`, the pipeline automatically runs the downstream R script:

```bash
Rscript "${SCRIPTS_DIR}/run_DNAcopy_from_bash.R" "path/to/sample1_sample2.ratio_per_w_CC0.log2adj.txt"
```

Alternatively, you can run this R step manually in your terminal:

```bash
input_for_r="path/to/sample1_sample2.ratio_per_w_CC0.log2adj.txt"
Rscript "${SCRIPTS_DIR}/run_DNAcopy_from_bash.R" "$input_for_r"
```

---

## Outputs

Results are organized under `OUTPUT_DIR`:

* `samtools_stats/` – Individual BAM stats files
* `multiqc/` – MultiQC HTML reports
* `merged_bams/` – Merged BAM files and indices
* `cov_stats/` – Modal and min/max coverage tables
* `difcover/` – Differential coverage ratios
* `dnacopy/` – log2-adjusted ratios for DNAcopy segmentation

---

## Notes

* Paired BAMs for each species (`Male`, `Female`) are expected.
* Paths in `config.yaml` must be customized for your environment.
* Logs for each rule are stored under `OUTPUT_DIR/logs/`.
* Ensure `samtools`, `multiqc`, and R (`DNAcopy`) are installed in your Conda environments.

---
