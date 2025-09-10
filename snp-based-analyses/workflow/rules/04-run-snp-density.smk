# Extract configuration parameters
import os

bcf_dir = config["INDV_BCF_DIR"]
threads = int(config["THREADS"])
output_dir = config["OUTPUT_DIR"]
species_list = config["SPECIES"]
if isinstance(species_list, str):
    species_list = [species_list]

# Define the output directory for SNP density
output_snp_dir = os.path.join(output_dir, "snpdensity")

# Function to collect sample names from BCF directory structure
def get_samples(species_list):
    samples = { sp: {"Male": [], "Female": []} for sp in species_list }
    for sp in species_list:
        sp_dir = os.path.join(bcf_dir, sp)
        for gender in ["Male", "Female"]:
            gender_dir = os.path.join(sp_dir, gender)
            if os.path.isdir(gender_dir):
                for f in os.listdir(gender_dir):
                    if f.endswith(".bcf"):
                        sample = os.path.splitext(f)[0]
                        samples[sp][gender].append(sample)
    return samples

# Get dictionary of all samples by species and gender
samples_dict = get_samples(species_list)

rule calculate_snp_density_male:
    input:
        bcf = lambda wc: f"{bcf_dir}/{wc.species}/Male/{wc.sample}.bcf"
    output:
        f"{output_snp_dir}/{{species}}/Male_{{sample}}_10k.snpden"
    params:
        out_prefix = lambda wc: f"{output_snp_dir}/{wc.species}/Male_{wc.sample}_10k"
    threads: threads
    shell:
        """
        module load vcftools/0.1.16
        vcftools --bcf {input.bcf} --SNPdensity 10000 --out {params.out_prefix}
        """

rule calculate_snp_density_female:
    input:
        bcf = lambda wc: f"{bcf_dir}/{wc.species}/Female/{wc.sample}.bcf"
    output:
        f"{output_snp_dir}/{{species}}/Female_{{sample}}_10k.snpden"
    params:
        out_prefix = lambda wc: f"{output_snp_dir}/{wc.species}/Female_{wc.sample}_10k"
    threads: threads
    shell:
        """
        module load vcftools/0.1.16
        vcftools --bcf {input.bcf} --SNPdensity 10000 --out {params.out_prefix}
        """
