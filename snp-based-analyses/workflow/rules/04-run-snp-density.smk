import os

# ─────────────────────────────────────────────────────────────
# Rule file: 04-run-snp-density.smk
# Calculate SNP density for individual sample BCF files
# ─────────────────────────────────────────────────────────────

# Configuration parameters
bcf_dir = config["INDV_BCF_DIR"]
output_dir = config["OUTPUT_DIR"]
threads = int(config.get("THREADS", 4))

# Define SNP density output base directory
output_snp_dir = os.path.join(output_dir, "snpdensity")

# ─────────────────────────────────────────────────────────────
# Helper: Get BCF input file for a given sample
# ─────────────────────────────────────────────────────────────
def get_bcf_for_sample(wildcards):
    """
    Return the BCF file path for a given species/gender/sample combination.
    This BCF was created by rule extract_individual_bcf.
    """
    return os.path.join(
        config["INDV_BCF_DIR"],
        wildcards.species,
        wildcards.gender,
        f"{wildcards.sample}.bcf"
    )

# ─────────────────────────────────────────────────────────────
# Rule: Calculate SNP density for each sample
# ─────────────────────────────────────────────────────────────
rule snp_density:
    input:
        bcf = get_bcf_for_sample,
        index = lambda wildcards: get_bcf_for_sample(wildcards) + ".csi"
    output:
        snpden = os.path.join(output_snp_dir, "{species}", "{gender,Male|Female}_{sample}_10k.snpden")
    params:
        out_prefix = lambda wildcards, output: output.snpden[:-7]
    threads: threads
    conda:
        "envs/fst_env.yaml"
    log:
        os.path.join(output_snp_dir, "{species}", "logs", "{gender}_{sample}_snpdensity.log")
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {output.snpden})
        mkdir -p $(dirname {log})
        
        echo "=== Starting SNP density calculation ===" > {log} 2>&1
        echo "Species: {wildcards.species}" >> {log} 2>&1
        echo "Gender: {wildcards.gender}" >> {log} 2>&1
        echo "Sample: {wildcards.sample}" >> {log} 2>&1
        echo "Input BCF: {input.bcf}" >> {log} 2>&1
        echo "Output: {output.snpden}" >> {log} 2>&1
        echo "========================================" >> {log} 2>&1
        
        # Verify input files exist
        if [ ! -f {input.bcf} ]; then
            echo "ERROR: Input BCF file not found: {input.bcf}" >> {log} 2>&1
            exit 1
        fi
        
        if [ ! -f {input.index} ]; then
            echo "ERROR: Input index file not found: {input.index}" >> {log} 2>&1
            exit 1
        fi
        
        # Verify vcftools is available
        echo "=== Verifying vcftools installation ===" >> {log} 2>&1
        which vcftools >> {log} 2>&1
        vcftools --version >> {log} 2>&1
        echo "========================================" >> {log} 2>&1
            
        # Run vcftools SNP density calculation
        echo "=== Running vcftools ===" >> {log} 2>&1
        vcftools --bcf {input.bcf} \
            --SNPdensity 10000 \
            --out {params.out_prefix} >> {log} 2>&1
        
        # Check vcftools exit status
        VCFTOOLS_EXIT=$?
        if [ $VCFTOOLS_EXIT -ne 0 ]; then
            echo "ERROR: vcftools failed with exit code $VCFTOOLS_EXIT" >> {log} 2>&1
            exit 1
        fi
        
        # Verify output
        if [ -f {output.snpden} ]; then
            echo "=== SUCCESS: Output file created ===" >> {log} 2>&1
            ls -lh {output.snpden} >> {log} 2>&1
            wc -l {output.snpden} >> {log} 2>&1
        else
            echo "=== ERROR: Output file not created ===" >> {log} 2>&1
            echo "Expected output: {output.snpden}" >> {log} 2>&1
            echo "Directory contents:" >> {log} 2>&1
            ls -lh $(dirname {output.snpden}) >> {log} 2>&1
            exit 1
        fi
        
        echo "=== SNP density calculation complete ===" >> {log} 2>&1
        """

# ─────────────────────────────────────────────────────────────
# Rule: Plot SNP density sex differences
# ─────────────────────────────────────────────────────────────
rule plot_snp_density:
    input:
        male_snpden=lambda wc: expand(
            os.path.join(output_snp_dir, "{species}", "Male_{sample}_10k.snpden"),
            species=[wc.species],
            sample=samples_dict[wc.species]["Male"]
        ),
        female_snpden=lambda wc: expand(
            os.path.join(output_snp_dir, "{species}", "Female_{sample}_10k.snpden"),
            species=[wc.species],
            sample=samples_dict[wc.species]["Female"]
        ),
        ref=config["ref_file"]
    output:
        plot=os.path.join(output_snp_dir, "{species}", "plot", "{species}_proportion_p001.pdf"),
        table=os.path.join(output_snp_dir, "{species}", "plot", "{species}_SNPdensity_SexFindR.txt"),
        sig_table=os.path.join(output_snp_dir, "{species}", "plot", "{species}_significant_SNPdensity_SexFindR.txt")
    params:
        script=config["snp_density_rscript"],
        snp_base_dir=output_snp_dir
    conda:
        "envs/r-plots_env.yaml"
    shell:
        """
        Rscript {params.script} {params.snp_base_dir} {input.ref} {wildcards.species}
        """