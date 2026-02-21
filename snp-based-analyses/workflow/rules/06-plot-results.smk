import os

# ─────────────────────────────────────────────────────────────
# Integrated Plotting and Analysis
# ─────────────────────────────────────────────────────────────
# This module creates integrated plots combining results from
# all four analyses: FST, nucleotide diversity, SNP density, and GWAS.
# It identifies candidate sex-linked regions that rank highly across
# multiple methods and calculates correlation between approaches.
# ─────────────────────────────────────────────────────────────

OUTPUT_DIR = config["OUTPUT_DIR"]

# ─────────────────────────────────────────────────────────────
# Rule: Integrated analysis and combined_plots
# ─────────────────────────────────────────────────────────────
rule integrated_analysis:
    input:
        snp_density=os.path.join(OUTPUT_DIR, "snpdensity", "{species}", "plot", "{species}_SNPdensity_SexFindR.txt"),
        gwas_window=os.path.join(OUTPUT_DIR, "snp_gwas", "gwas", "{species}", "output", "{species}_GWAS_window_count.txt"),
        fst_window=os.path.join(OUTPUT_DIR, "fst", "{species}", "{species}_Fst_window_count.txt"),
        pi_diff=os.path.join(OUTPUT_DIR, "pi", "{species}", "{species}_NucDiv_10kb_sex_differences.txt"),
        ref=config["ref_file"]
    output:
        candidates_all=os.path.join(OUTPUT_DIR, "combined_plots", "{species}", "candidate_windows_SNP_GWAS_Fst_pi.txt"),
        candidates_gwas_pi_fst=os.path.join(OUTPUT_DIR, "combined_plots", "{species}", "candidate_windows_GWAS_Pi_Fst.txt"),
        candidates_snp_fst_pi=os.path.join(OUTPUT_DIR, "combined_plots", "{species}", "candidate_windows_SNP_Fst_Pi.txt"),
        correlations=os.path.join(OUTPUT_DIR, "combined_plots", "{species}", "{species}_correlation_rankings.txt"),
        plot=os.path.join(OUTPUT_DIR, "combined_plots", "{species}", "{species}_integrated_plot.png")
    params:
        script=config["integrated_analysis_rscript"]
    conda:
        "envs/r-plots_env.yaml"
    shell:
        """
        mkdir -p $(dirname {output.plot})
        Rscript {params.script} {input.snp_density} {input.gwas_window} {input.fst_window} {input.pi_diff} {input.ref} {wildcards.species} {output.candidates_all} {output.candidates_gwas_pi_fst} {output.candidates_snp_fst_pi} {output.correlations} {output.plot}
        """