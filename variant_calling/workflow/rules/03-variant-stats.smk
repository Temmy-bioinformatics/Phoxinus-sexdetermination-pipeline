########################################
########## Load Configuration ##########
########################################

configfile: "config/config.yaml"  # Replace with the actual filename

SAMPLE = config["SAMPLE_ID"]
BASE_OUTPUT_DIR = config["BASE_OUTPUT_DIR"].format(SAMPLE_ID=SAMPLE)
PLOT_DIR = f"{BASE_OUTPUT_DIR}/Plots"
SUBSET_VCF = config["DERIVED_PATHS"]["MAF_FINAL_OUTPUT"] + ".vcf"
OUT_PREFIX = f"{BASE_OUTPUT_DIR}/Stats/{SAMPLE}"

########################################
############## Rules ###################
########################################

rule all:
    input:
        f"{OUT_PREFIX}.frq",
        f"{OUT_PREFIX}.idepth",
        f"{OUT_PREFIX}.ldepth.mean",
        f"{OUT_PREFIX}.lqual",
        f"{OUT_PREFIX}.imiss",
        f"{OUT_PREFIX}.lmiss",
        f"{OUT_PREFIX}.het",
        f"{PLOT_DIR}/variants_statistics.png"

rule calculate_allele_frequency:
    input:
        vcf=SUBSET_VCF
    output:
        freq=f"{OUT_PREFIX}.frq"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input.vcf} --freq2 --out {OUT_PREFIX}"

rule calculate_mean_depth_individual:
    input:
        vcf=SUBSET_VCF
    output:
        depth_indv=f"{OUT_PREFIX}.idepth"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input.vcf} --depth --out {OUT_PREFIX}"

rule calculate_mean_depth_site:
    input:
        vcf=SUBSET_VCF
    output:
        depth_site=f"{OUT_PREFIX}.ldepth.mean"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input.vcf} --site-mean-depth --out {OUT_PREFIX}"

rule calculate_site_quality:
    input:
        vcf=SUBSET_VCF
    output:
        site_quality=f"{OUT_PREFIX}.lqual"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input.vcf} --site-quality --out {OUT_PREFIX}"

rule calculate_missing_data_individual:
    input:
        vcf=SUBSET_VCF
    output:
        missing_indv=f"{OUT_PREFIX}.imiss"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input.vcf} --missing-indv --out {OUT_PREFIX}"

rule calculate_missing_data_site:
    input:
        vcf=SUBSET_VCF
    output:
        missing_site=f"{OUT_PREFIX}.lmiss"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input.vcf} --missing-site --out {OUT_PREFIX}"

rule calculate_f_statistics:
    input:
        vcf=SUBSET_VCF
    output:
        het=f"{OUT_PREFIX}.het"
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --vcf {input.vcf} --het --out {OUT_PREFIX}"

rule generate_plots:
    input:
        frq=f"{OUT_PREFIX}.frq",
        indv=f"{OUT_PREFIX}.idepth",
        ldepth_mean=f"{OUT_PREFIX}.ldepth.mean",
        lqual=f"{OUT_PREFIX}.lqual",
        imiss=f"{OUT_PREFIX}.imiss",
        lmiss=f"{OUT_PREFIX}.lmiss",
        het=f"{OUT_PREFIX}.het"
    output:
        f"{PLOT_DIR}/variants_statistics.png"
    conda:
        "envs/R_env.yml"
    script:
        "workflow/rscripts/Plot_variant_stats.r"
