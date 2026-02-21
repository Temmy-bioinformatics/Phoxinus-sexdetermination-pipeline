import os

# ─────────────────────────────────────────────────────────────
# FST Analysis Rules
# ─────────────────────────────────────────────────────────────
# This module calculates FST (fixation index) between male and 
# female populations using VCFtools, both per-SNP and in 10kb windows.
# It also generates Manhattan-style plots of windowed FST values.
# ─────────────────────────────────────────────────────────────

# Use paths from config
OUTPUT_DIR = config["OUTPUT_DIR"]
BCF_DIR = config["BCF_DIR"]
SAMPLE_DIR = config["SAMPLE_DIR"]

# ─────────────────────────────────────────────────────────────
# Rule: Calculate per-SNP FST
# ─────────────────────────────────────────────────────────────
rule run_fst:
    input:
        bcf=lambda wc: f"{BCF_DIR}/{wc.species}/Filtered/{wc.species}_allChromosomes.filtered.norm.bcf",
        pop1=lambda wc: f"{SAMPLE_DIR}/{wc.species}_Male.vcflist",
        pop2=lambda wc: f"{SAMPLE_DIR}/{wc.species}_Female.vcflist"
    output:
        fst=os.path.join(
            OUTPUT_DIR, "fst", "{species}", "{species}_Fst.txt.weir.fst"
        )
    params:
        out_prefix=os.path.join(
            OUTPUT_DIR, "fst", "{species}", "{species}_Fst.txt"
        )
    conda:
        "envs/fst_env.yaml"
    shell:
        """
        mkdir -p $(dirname {output.fst})
        vcftools --bcf {input.bcf} \
                 --weir-fst-pop {input.pop1} \
                 --weir-fst-pop {input.pop2} \
                 --out {params.out_prefix}
        """

# ─────────────────────────────────────────────────────────────
# Rule: Calculate FST in 10kb windows
# ─────────────────────────────────────────────────────────────
rule run_fst_10k_windows:
    input:
        bcf=lambda wc: f"{BCF_DIR}/{wc.species}/Filtered/{wc.species}_allChromosomes.filtered.norm.bcf",
        pop1=lambda wc: f"{SAMPLE_DIR}/{wc.species}_Male.vcflist",
        pop2=lambda wc: f"{SAMPLE_DIR}/{wc.species}_Female.vcflist"
    output:
        fst=os.path.join(
            OUTPUT_DIR, "fst", "{species}", "{species}_Fst_10k.windowed.weir.fst"
        )
    params:
        out_prefix=os.path.join(
            OUTPUT_DIR, "fst", "{species}", "{species}_Fst_10k"
        )
    conda:
        "envs/fst_env.yaml"
    shell:
        """
        mkdir -p $(dirname {output.fst})
        vcftools --bcf {input.bcf} \
                 --weir-fst-pop {input.pop1} \
                 --weir-fst-pop {input.pop2} \
                 --fst-window-size 10000 \
                 --out {params.out_prefix}
        """

# ─────────────────────────────────────────────────────────────
# Rule: Extract top 5% FST SNPs from per-SNP data
# ─────────────────────────────────────────────────────────────
rule extract_top_fst_snps:
    input:
        fst=os.path.join(
            OUTPUT_DIR, "fst", "{species}", "{species}_Fst.txt.weir.fst"
        ),
        ref=config["ref_file"]
    output:
        snps=os.path.join(
            OUTPUT_DIR, "fst", "{species}", "{species}_all_Fst_sort_cut.txt"
        )
    params:
        script=config["extract_fst_rscript"]
    conda:
        "envs/r-plots_env.yaml"
    shell:
        """
        Rscript {params.script} {input.fst} {input.ref} {wildcards.species} {output.snps}
        """

# ─────────────────────────────────────────────────────────────
# Rule: Sort extracted FST SNPs by coordinates
# ─────────────────────────────────────────────────────────────
rule sort_top_fst_snps:
    input:
        snps=os.path.join(
            OUTPUT_DIR, "fst", "{species}", "{species}_all_Fst_sort_cut.txt"
        )
    output:
        sorted=os.path.join(
            OUTPUT_DIR, "fst", "{species}", "{species}_coordinated_sorted_Fst_cut.txt"
        )
    shell:
        """
        sort -k1,1 -k2n {input.snps} > {output.sorted}
        """

# ─────────────────────────────────────────────────────────────
# Rule: Count top FST SNPs per 10kb window
# ─────────────────────────────────────────────────────────────
rule count_fst_windows:
    input:
        sorted=os.path.join(
            OUTPUT_DIR, "fst", "{species}", "{species}_coordinated_sorted_Fst_cut.txt"
        )
    output:
        counts=os.path.join(
            OUTPUT_DIR, "fst", "{species}", "{species}_Fst_window_count.txt"
        )
    params:
        script=config["count_window"],
        window_size=10000
    shell:
        """
        python {params.script} {input.sorted} {output.counts} {params.window_size}
        """

# ─────────────────────────────────────────────────────────────
# Rule: Plot 10kb windowed FST
# ─────────────────────────────────────────────────────────────
rule plot_fst:
    input:
        fst=os.path.join(
            OUTPUT_DIR, "fst", "{species}", "{species}_Fst_10k.windowed.weir.fst"
        ),
        ref=config["ref_file"]
    output:
        plot=os.path.join(
            OUTPUT_DIR, "plotting", "{species}", "{species}_Fst_10k_plot.tiff"
        ),
        thresholds=os.path.join(
            OUTPUT_DIR, "plotting", "{species}", "{species}_Fst_10k_thresholds.txt"
        )
    params:
        script=config["fst_rscript"]
    conda:
        "envs/r-plots_env.yaml"
    shell:
        """
        mkdir -p $(dirname {output.plot})
        Rscript {params.script} {input.fst} {input.ref} {wildcards.species} {output.plot} {output.thresholds}
        """
