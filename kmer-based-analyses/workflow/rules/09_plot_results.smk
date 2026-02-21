# ─────────────────────────────────────────────────────────────
# K-mer BLAST Plotting Rule
# ─────────────────────────────────────────────────────────────
# This module creates chromosome-wide visualizations of k-mer
# BLAST results from module 07 (RUN_BLAST_KMERS).
# ─────────────────────────────────────────────────────────────

import os

PROJECT_DIR = config["PROJECT_DIR"]
SPECIES = config["SPECIES"]

# ─────────────────────────────────────────────────────────────
# Rule: Plot k-mer BLAST results (chromosome-wide only)
# ─────────────────────────────────────────────────────────────
rule plot_kmer_blast:
    input:
        blast=os.path.join(PROJECT_DIR, SPECIES, "blast", "kmers_blast.out"),
        ref=config["REFERENCE"] + ".fai"
    output:
        chr_plot_png=os.path.join(
            PROJECT_DIR, SPECIES, "plots",
            SPECIES + "_kmer_counts_per_chr.png"
        ),
        chr_plot_pdf=os.path.join(
            PROJECT_DIR, SPECIES, "plots",
            SPECIES + "_kmer_counts_per_chr.pdf"
        ),
        chr_plot_svg=os.path.join(
            PROJECT_DIR, SPECIES, "plots",
            SPECIES + "_kmer_counts_per_chr.svg"
        ),
        summary=os.path.join(
            PROJECT_DIR, SPECIES, "plots",
            SPECIES + "_kmer_blast_mapping_summary.csv"
        )
    params:
        script=config["plot_kmer_blast_script"],
        species_name=config.get("SPECIES_NAMES", {}).get(SPECIES, SPECIES),
        output_dir=os.path.join(PROJECT_DIR, SPECIES, "plots")
    conda:
        "envs/plot_blast_results.yaml"
    shell:
        """
        Rscript {params.script} \
            {input.blast} \
            {input.ref} \
            "{params.species_name}" \
            {params.output_dir} \
            {SPECIES}
        """
