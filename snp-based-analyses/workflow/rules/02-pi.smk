# Extract variables from config
project_dir = config["PROJECT_DIR"]
bcf_dir = config["BCF_DIR"]
sample_dir = config["SAMPLE_DIR"]
output_base_dir = config["OUTPUT_DIR"]
vcf_lists = config["VCF_LISTS"]
threads = config["THREADS"]

# Extract species names from the list of vcflist files
species_files = [f.split("_")[0] for f in vcf_lists]
species = sorted(list(set(species_files)))
genders = ["Male", "Female"]


rule run_nucleotide_diversity:
    input:
        bcf=lambda wildcards: f"{bcf_dir}/{wildcards.species}/Filtered/{wildcards.species}_allChromosomes.filtered.norm.bcf",
        vcf_list=lambda wildcards: f"{sample_dir}/{wildcards.species}_{wildcards.gender}.vcflist"
    output:
        diversity = "{output_base_dir}/pi/{species}/{species}_{gender}_10k_nucleotide_diversity.windowed.pi"
    threads: threads
    #conda:
    #    "./envs/fst_env.yaml"
    shell:
        """
        mkdir -p $(dirname {output.diversity})
        OUT_PREFIX=$(echo {output.diversity} | sed 's/.windowed.pi//')
        vcftools --bcf {input.bcf} \
                 --keep {input.vcf_list} \
                 --window-pi 10000 \
                 --out $OUT_PREFIX
        """
