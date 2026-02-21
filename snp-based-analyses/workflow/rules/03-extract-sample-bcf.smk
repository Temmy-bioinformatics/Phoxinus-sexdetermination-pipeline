# ─────────────────────────────────────────────────────────────
# Rule file: 03-extract-sample-bcf.smk
# Extract individual sample BCF files from merged BCF
# ─────────────────────────────────────────────────────────────

# Configuration parameters
bcf_dir = config["BCF_DIR"]
indv_bcf_dir = config["INDV_BCF_DIR"]
threads = config.get("THREADS", 4)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
# ─────────────────────────────────────────────────────────────
# Rule: Extract individual BCF files for all samples
# ─────────────────────────────────────────────────────────────
rule extract_individual_bcf:
    input:
        bcf = lambda wildcards: os.path.join(
            bcf_dir,
            wildcards.species,
            "Filtered",
            f"{wildcards.species}_allChromosomes.filtered.norm.bcf"
        )
    output:
        bcf = os.path.join(indv_bcf_dir, "{species}", "{gender}", "{sample}.bcf"),
        index = os.path.join(indv_bcf_dir, "{species}", "{gender}", "{sample}.bcf.csi")
    threads: threads
    conda:
        "/home/toriowo/PROJECTS/sex-determination-pipeline/code/job-scripts/snp-based-analyses/workflow/envs/bcftools_env.yaml"
    shell:
        """
        # Verify bcftools installation
        echo "=== Verifying bcftools installation ==="
        which bcftools
        bcftools --version
        echo "======================================="
        
        mkdir -p $(dirname {output.bcf})
        bcftools view -a -s {wildcards.sample} --threads {threads} -O b -o {output.bcf} {input.bcf}
        bcftools index {output.bcf}
        
        # Verify output was created successfully
        echo "=== Verifying output ==="
        ls -lh {output.bcf}
        ls -lh {output.index}
        echo "========================"
        """
