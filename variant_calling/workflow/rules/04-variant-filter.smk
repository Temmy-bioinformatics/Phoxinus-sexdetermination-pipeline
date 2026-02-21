########################################
########## Load Configuration ##########
########################################

configfile: "config/config.yaml" 

SAMPLE = config["SAMPLE_ID"]
THREADS = config["THREADS"]

# Paths from config
INPUT_BCF = config["DERIVED_PATHS"]["FILTERED_BCF"]
NORM_BCF = config["DERIVED_PATHS"]["NORM_BCF"]

# Filters
QUAL = config["QUAL"]
MIN_DEPTH = config["MIN_DEPTH"]
MAX_DEPTH = config["MAX_DEPTH"]
F_MISSING = config["F_MISSING"]

########################################
############## Rules ###################
########################################

rule all:
    input:
        NORM_BCF

rule bcftools_filter:
    input:
        bcf=INPUT_BCF
    output:
        filtered_bcf=NORM_BCF
    threads: THREADS
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        bcftools view -i '(QUAL>={QUAL} & FORMAT/DP>={MIN_DEPTH} & FORMAT/DP<={MAX_DEPTH}) & TYPE!="indel"' \
            --threads {threads} -Ou {input.bcf} | \
        bcftools filter -i 'F_MISSING<={F_MISSING}' --threads {threads} -Ou | \
        bcftools norm -d all --threads {threads} --write-index -Ob -o {output.filtered_bcf}
        """
