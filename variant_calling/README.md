# Snakemake Pipeline Configuration

This Snakemake pipeline is designed for variant calling, subsampling, statistics, and filtering steps in genomic analyses. 

## Configuration and Setup

- **Config file:** The pipeline uses a YAML config file (`config/{your_config_filename}.yaml`) to specify sample IDs, threads, output directories, and pipeline parameters.
- **Chromosomes:** By default, chromosomes `Chr1` to `Chr25` are included but can be adjusted.
- **Toggle switches:** The pipeline steps can be enabled or disabled using these boolean flags in the config:
  - `RUN_CALLING` – Variant calling with bcftools
  - `RUN_SUBSAMPLE` – Subsampling variants
  - `RUN_STATS` – Generate variant statistics and plots
  - `RUN_FILTER` – Variant filtering

## Included Modules

The pipeline conditionally includes rules based on the toggles:

- `01-bcftools-call-variants.smk`
- `02-subsample-variants.smk`
- `03-variant-stats.smk`
- `04-variant-filter.smk`

## Final Targets

The `rule all` aggregates output files based on enabled modules, including:

- Merged chromosome BCF files
- Subsampled VCF file
- Variant statistics files and plots
- Filtered BCF file

## Usage

1. Update your config YAML file with sample-specific parameters.
2. Run Snakemake with this master Snakefile.
3. Outputs will be organized in directories under the base output path.

---

This setup enables flexible and modular variant processing for scalable genomic workflows.
