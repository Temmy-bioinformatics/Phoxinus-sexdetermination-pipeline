##############################
######### Snakemake ##########
##############################

# Load configuration file
configfile: "config/config.yaml"  # Update with actual config name

# Extract variables
SAMPLE = config["SAMPLE_ID"]
THREADS = config["THREADS"]
MAXMEMORY = config["MAXMEMORY"]
BASE_OUTPUT_DIR = config["BASE_OUTPUT_DIR"].format(SAMPLE_ID=SAMPLE)
SUBSAMPLE_SIZE = config["SUBSAMPLE_SIZE"]

# Paths
MERGED_BCF = f"{BASE_OUTPUT_DIR}/{SAMPLE}_allChromosomes_merged.bcf"
SUBSAMPLED_VCF = f"{BASE_OUTPUT_DIR}/Subsample_vcf/Subsampled_{SUBSAMPLE_SIZE}.vcf"
HEADER_FILE = f"{BASE_OUTPUT_DIR}/Subsample_vcf/header_{SUBSAMPLE_SIZE}.txt"
VARIANTS_FILE = f"{BASE_OUTPUT_DIR}/Subsample_vcf/Subsampled_{SUBSAMPLE_SIZE}_Delimitation_variants.txt"

##############################
######### Rule Section #######
##############################

rule all:
    input:
        SUBSAMPLED_VCF

rule subsample_variants:
    input:
        bcf=MERGED_BCF,
        bcf_index=f"{MERGED_BCF}.csi"
    output:
        vcf=SUBSAMPLED_VCF
    params:
        header_file=HEADER_FILE,
        variants_file=VARIANTS_FILE,
        subsample_size=SUBSAMPLE_SIZE
    threads: THREADS
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        mkdir -p $(dirname {output.vcf})

        # Extract the header
        bcftools view -h {input.bcf} > {params.header_file}

        # Subsample N variants randomly
        bcftools view -H {input.bcf} --threads {threads} | shuf -n {params.subsample_size} > {params.variants_file}

        # Combine header and variants
        cat {params.header_file} {params.variants_file} > {output.vcf}
        """
