# Extract configuration parameters
project_dir = config["PROJECT_DIR"]
bcf_dir = config["BCF_DIR"]
sample_dir = config["SAMPLE_DIR"]
vcf_lists = config["VCF_LISTS"]
threads = config["THREADS"]
output_dir = config["OUTPUT_DIR"]

# Define subdirectory for output BCF files
output_vcf_dir = f"{output_dir}/individual_bcf"

# Parse species and genders from the VCF lists
species_files = [f.split("_")[0] for f in vcf_lists]
species = sorted(set(species_files))

# Function to parse samples from VCF list files
def get_samples(vcf_lists, gender):
    samples = {}
    for vcf_list in vcf_lists:
        sp, vcf_gender = vcf_list.rsplit("_", 1)
        vcf_gender = vcf_gender.split(".")[0]
        if vcf_gender == gender:
            with open(f"{sample_dir}/{vcf_list}", 'r') as file:
                samples[sp] = file.read().splitlines()
    return samples

# Retrieve sample dicts
male_samples_dict = get_samples(vcf_lists, "Male")
female_samples_dict = get_samples(vcf_lists, "Female")

rule extract_individual_bcfs_male:
    input:
        bcf = lambda wildcards: f"{bcf_dir}/{wildcards.species}/Filtered/{wildcards.species}_allChromosomes.filtered.norm.bcf"
    output:
        f"{output_vcf_dir}/{{species}}/Male/{{sample}}.bcf"
    threads: threads
    shell:
        """
        mkdir -p $(dirname {output})
        bcftools view -a -s {wildcards.sample} --threads {threads} -O b -o {output} {input.bcf}
        bcftools index {output}
        """

rule extract_individual_bcfs_female:
    input:
        bcf = lambda wildcards: f"{bcf_dir}/{wildcards.species}/Filtered/{wildcards.species}_allChromosomes.filtered.norm.bcf"
    output:
        f"{output_vcf_dir}/{{species}}/Female/{{sample}}.bcf"
    threads: threads
    shell:
        """
        mkdir -p $(dirname {output})
        bcftools view -a -s {wildcards.sample} --threads {threads} -O b -o {output} {input.bcf}
        bcftools index {output}
        """
