##############################
######### Snakemake ##########
##############################

# Manually define number of chromosomes
Chrs = [f"Chr{i}" for i in range(1, 26)]

# Load config
configfile: "config/config.yaml"

# Extract configuration variables
SAMPLE = config["SAMPLE_ID"]
REFGENOME = config["REFGENOME"]
BAM_LIST = config["BAM_LIST"].format(SAMPLE_ID=SAMPLE)
BASE_OUTPUT_DIR = config["BASE_OUTPUT_DIR"].format(SAMPLE_ID=SAMPLE)
THREADS = config["THREADS"]
MAXMEMORY = config["MAXMEMORY"]

# Derived paths
SPLIT_DIR = f"{BASE_OUTPUT_DIR}/Split_Chrs"
LOG_DIR = f"{BASE_OUTPUT_DIR}/logs"
SUBSAMPLED_VCF = f"{BASE_OUTPUT_DIR}/Subsample_vcf/Subsampled_10k.vcf"
MERGED_BCF = f"{BASE_OUTPUT_DIR}/{SAMPLE}_allChromosomes_merged.bcf"

##############################
###### Rules Section #########
##############################

rule all:
    input:
        expand(f"{SPLIT_DIR}/{SAMPLE}_{{chr}}_merged.bcf", chr=Chrs),
        MERGED_BCF,
        SUBSAMPLED_VCF

rule bcftools_call:
    input:
        ref=REFGENOME,
        bam_list=BAM_LIST
    output:
        bcf=f"{SPLIT_DIR}/{SAMPLE}_{{chr}}_merged.bcf",
        bcf_index=f"{SPLIT_DIR}/{SAMPLE}_{{chr}}_merged.bcf.csi"
    params:
        chr=lambda wildcards: wildcards.chr,
        start_log=f"{LOG_DIR}/{{chr}}_start.log",
        end_log=f"{LOG_DIR}/{{chr}}_end.log"
    threads: THREADS
    resources:
        mem=MAXMEMORY
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        bcftools mpileup -q 20 -Q 20 -Ou \
          -a FORMAT/DP,AD,INFO/AD -P ILLUMINA -f {input.ref} \
          -b {input.bam_list} -r {params.chr} | \
        bcftools call -f GQ -mv -Ob -o {output.bcf}

        bcftools index {output.bcf}
        """

rule bcftools_concat:
    input:
        expand(f"{SPLIT_DIR}/{SAMPLE}_{{chr}}_merged.bcf", chr=Chrs)
    output:
        bcf=MERGED_BCF,
        bcf_index=f"{MERGED_BCF}.csi"
    threads: THREADS
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        bcftools concat {input} --threads {threads} --write-index -Ob --output {output.bcf}
        """
