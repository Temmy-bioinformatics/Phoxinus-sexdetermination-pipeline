import os

# Get species list by finding all *_Female.vcflist files
species_files = [f for f in os.listdir(config["SAMPLE_DIR"]) if f.endswith("_Female.vcflist")]
species = sorted([f.replace("_Female.vcflist", "") for f in species_files])

# Rule: Regular Fst calculation
rule run_fst:
    input:
        bcf=lambda wc: f"{config['BCF_DIR']}/{wc.species}/Filtered/{wc.species}_allChromosomes.filtered.norm.bcf",
        pop1=lambda wc: f"{config['SAMPLE_DIR']}/{wc.species}_Male.vcflist",
        pop2=lambda wc: f"{config['SAMPLE_DIR']}/{wc.species}_Female.vcflist"
    output:
        fst="results/fst/{species}/{species}_Fst.txt.weir.fst"
    params:
        out_prefix="results/fst/{species}/{species}_Fst.txt"
    #conda:
    #    "/home/toriowo/Snakemake_projects/sexfindr-snakemake/snp-based-analyses/workflow/envs/fst_env.yaml"
    shell:
        """
        mkdir -p $(dirname {output.fst})
        vcftools --bcf {input.bcf} \
                 --weir-fst-pop {input.pop1} \
                 --weir-fst-pop {input.pop2} \
                 --out {params.out_prefix}
        """

# Rule: 10k windowed Fst calculation
rule run_fst_10k_windows:
    input:
        bcf=lambda wc: f"{config['BCF_DIR']}/{wc.species}/Filtered/{wc.species}_allChromosomes.filtered.norm.bcf",
        pop1=lambda wc: f"{config['SAMPLE_DIR']}/{wc.species}_Male.vcflist",
        pop2=lambda wc: f"{config['SAMPLE_DIR']}/{wc.species}_Female.vcflist"
    output:
        fst="results/fst/{species}/{species}_Fst_10k.windowed.weir.fst"
    params:
        out_prefix="results/fst/{species}/{species}_Fst_10k"
    #conda:
    #    "/home/toriowo/Snakemake_projects/sexfindr-snakemake/snp-based-analyses/workflow/envs/fst_env.yaml"
    shell:
        """
        mkdir -p $(dirname {output.fst})
        vcftools --bcf {input.bcf} \
                 --weir-fst-pop {input.pop1} \
                 --weir-fst-pop {input.pop2} \
                 --fst-window-size 10000 \
                 --out {params.out_prefix}
        """
