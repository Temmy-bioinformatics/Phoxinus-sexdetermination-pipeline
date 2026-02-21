Coverage-Based Sex Determination Pipeline

This Snakemake pipeline performs coverage-based analyses for sex determination across multiple species. It includes BAM statistics, coverage analysis, differential coverage, and DNAcopy segmentation with plotting.

Pipeline Overview

The workflow consists of the following steps:

1. Samtools Stats (01-samtools-stats.smk)

Calculates BAM statistics for each species and sex using samtools stats.

Input: BAM metadata files listing all BAMs for a species and sex.

Output: Directory samtools_stats/{species}/{sex}/ containing .stats.txt for each BAM.

Logs: Stored in logs/samtools_stats/.

Environment: envs/samtools.yaml

2. MultiQC Report (02-multiqc.smk)

Aggregates all samtools stats into a single HTML report per species.

Input: samtools_stats/{species}/ (all sexes)

Output: multiqc/{species}_multiqc_report.html

Environment: envs/multiqc.yaml

3. Merge BAMs (03-merge-bams.smk)

Merges individual BAMs per species and sex.

Input: BAM list files from resources/bamlists/.

Output:

Merged BAM: {MERGED_BAM_DIR}/{species}/{species}_{sex}_merged.bam

BAM index: {MERGED_BAM_DIR}/{species}/{species}_{sex}_merged.bam.bai

Merged BAM stats: {MERGED_BAM_DIR}/{species}/samtools_stats/{species}_{sex}_merged.stats.txt

Environment: envs/samtools.yaml

Threads: 8

4. Coverage Stats (04-cov-stats.smk)

Computes min, max, mean coverage per sex.

Extracts modal coverage from samtools stats.

Calculates coverage ratio (male/female) and BAM size ratio.

Input: Merged BAMs and their stats.

Output: cov_stats/{species}/cov_stats_{species}.tab

Environment: envs/samtools.yaml

5. Differential Coverage (05-difcover.smk)

Computes per-window coverage ratios between sexes.

Uses coverage stats to set thresholds.

Runs:

from_bams_to_unionbed.sh → generates .unionbedcv

from_unionbed_to_ratio_per_window_CC0 → generates per-window ratios

Input:

Female and male merged BAMs

cov_stats_{species}.tab

Output:

difcover/{species}/sample1_sample2.unionbedcv

difcover/{species}/sample1_sample2.ratio_per_w_CC0.txt

Environment: envs/bedtools.yaml

6. DNACopy Input & Segmentation (06-dnacopy.smk + 07-dnacopy.smk)

Processes differential coverage ratios with log2 transformation using AC (coverage adjustment factor from cov_stats).

Runs DNAcopy segmentation in R.

Input: sample1_sample2.ratio_per_w_CC0.log2adj.txt

Output:

Segmentation results: dnacopy/{species}/sample1_sample2.ratio_per_w_CC0.log2adj.txt.DNAcopyout

Plots: dnacopy/{species}/sample1_sample2.ratio_per_w_CC0.log2adj.txt.pdf

Environment: envs/dna-copy-r.yaml

Script used: run_DNAcopy_from_bash.R in SCRIPTS_DIR

Configuration (config/config.yaml)

Key parameters:

PROJECT_DIR: Base project directory

OUTPUT_DIR: Results directory

MERGED_BAM_DIR: Directory to store merged BAMs

SPECIES: List of species to analyze

RUN_*: Toggles for each module (default: True)

SCRIPTS_DIR: Path to custom scripts

Running the Pipeline
snakemake --use-conda --cores <N>


Replace <N> with the number of CPU cores available.

The pipeline respects toggles in the configuration file.

Only enabled modules will run.

Directory Structure
PROJECT_DIR/
├── resources/bamlists/
│   └── {species}_{sex}_Bamfiles_Metadata.txt
├── results/
│   ├── samtools_stats/
│   ├── multiqc/
│   ├── merged_bams/
│   ├── cov_stats/
│   ├── difcover/
│   └── dnacopy/
├── workflow/
│   └── rules/
│       ├── 01-samtools-stats.smk
│       ├── 02-multiqc.smk
│       ├── 03-merge-bams.smk
│       ├── 04-cov-stats.smk
│       ├── 05-difcover.smk
│       ├── 06-dnacopy.smk
│       └── 07-dnacopy.smk
└── config/
    └── config.yaml

Outputs Summary
Step	Output files	Example paths
Samtools Stats	.stats.txt per BAM	results/samtools_stats/SeinePhoxinus/Female/sample.bam.stats.txt
MultiQC	HTML report	results/multiqc/SeinePhoxinus_multiqc_report.html
Merge BAMs	Merged BAM & index	merged_bams/SeinePhoxinus/SeinePhoxinus_Female_merged.bam
Coverage Stats	.tab summary	results/cov_stats/SeinePhoxinus/cov_stats_SeinePhoxinus.tab
Differential Coverage	.unionbedcv, .ratio_per_w_CC0.txt	results/difcover/SeinePhoxinus/sample1_sample2.ratio_per_w_CC0.txt
DNACopy	.DNAcopyout, .pdf	results/dnacopy/SeinePhoxinus/sample1_sample2.ratio_per_w_CC0.log2adj.txt.DNAcopyout
results/dnacopy/SeinePhoxinus/sample1_sample2.ratio_per_w_CC0.log2adj.txt.pdf
Notes

Ensure BAM files exist and are properly listed in metadata files.

Conda environments handle all software dependencies.

Directories are automatically created by the workflow.

Plots and output files are species-specific.

The pipeline is modular; disable steps using the RUN_* toggles in the config file.