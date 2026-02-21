import os

# ─────────────────────────────────────────────────────────────
# Config parameters
# ─────────────────────────────────────────────────────────────
OUTPUT_DIR = config["OUTPUT_DIR"]
BCF_DIR = config["BCF_DIR"]
SAMPLE_DIR = config["SAMPLE_DIR"]
THREADS = int(config.get("THREADS", 1))

# ─────────────────────────────────────────────────────────────
# Rule: Nucleotide diversity (Pi) calculation
# ─────────────────────────────────────────────────────────────
rule run_nucleotide_diversity:
    input:
        bcf=lambda wc: f"{BCF_DIR}/{wc.species}/Filtered/{wc.species}_allChromosomes.filtered.norm.bcf",
        vcf_list=lambda wc: f"{SAMPLE_DIR}/{wc.species}_{wc.gender}.vcflist"
    output:
        # Use a fixed string with wildcards; no lambda or dynamic
        "{output_dir}/pi/{species}/{species}_{gender}_10k_nucleotide_diversity.windowed.pi"
    threads: THREADS
    conda:
        "envs/fst_env.yaml"
    shell:
        """
        mkdir -p $(dirname {output})
        OUT_PREFIX=$(echo {output} | sed 's/.windowed.pi//')
        vcftools --bcf {input.bcf} \
                 --keep {input.vcf_list} \
                 --window-pi 10000 \
                 --out $OUT_PREFIX
        """


rule plot_nucleotide_diversity:
    input:
        male_pi=os.path.join(
            OUTPUT_DIR, "pi", "{species}", "{species}_Male_10k_nucleotide_diversity.windowed.pi"
        ),
        female_pi=os.path.join(
            OUTPUT_DIR, "pi", "{species}", "{species}_Female_10k_nucleotide_diversity.windowed.pi"
        ),
        ref=config["ref_file"]
    output:
        plot=os.path.join(
            OUTPUT_DIR, "pi", "{species}", "{species}_pi_10kb_window.tiff"
        ),
        var_plot=os.path.join(
            OUTPUT_DIR, "pi", "{species}", "{species}_var_per_10kb_window.tiff"
        ),
        diff_table=os.path.join(
            OUTPUT_DIR, "pi", "{species}", "{species}_NucDiv_10kb_sex_differences.txt"
        )
    params:
        script=config["nuc_diversity_rscript"],
        pi_dir=OUTPUT_DIR + "/pi"
    conda:
        "envs/r-plots_env.yaml"
    shell:
        """
        mkdir -p $(dirname {output.plot})
        Rscript {params.script} {params.pi_dir} {input.ref} {wildcards.species}
        """