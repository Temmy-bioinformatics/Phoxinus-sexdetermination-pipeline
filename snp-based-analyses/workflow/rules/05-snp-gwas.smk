import os

project_dir = config["PROJECT_DIR"]
bcf_dir = config["BCF_DIR"]
output_dir = config["OUTPUT_DIR"]
sample_list_dir = config["SAMPLE_DIR"]
threads = config["THREADS"]

bed_output_dir = os.path.join(output_dir, "snp_gwas", "bed_files")
gwas_output_dir = os.path.join(output_dir, "snp_gwas", "gwas")
pheno_output_dir = os.path.join(output_dir, "snp_gwas", "pheno_files")

# Rule to create phenotype file
rule create_pheno:
    input:
        male=lambda wc: f"{sample_list_dir}/{wc.species}_Male.vcflist",
        female=lambda wc: f"{sample_list_dir}/{wc.species}_Female.vcflist"
    output:
        pheno="{pheno_output_dir}/{species}/{species}.pheno.txt"
    threads: threads
    shell:
        """
        > {output.pheno}
        while read -r line; do
            echo -e "${{line}}\t${{line}}\t1" >> {output.pheno}
        done < {input.male}
        while read -r line; do
            echo -e "${{line}}\t${{line}}\t2" >> {output.pheno}
        done < {input.female}
        """

# Rule to create PLINK bed/bim/fam files
rule create_bed:
    input:
        bcf=lambda wc: f"{bcf_dir}/{wc.species}/Filtered/{wc.species}_allChromosomes.filtered.norm.bcf",
        pheno=lambda wc: f"{pheno_output_dir}/{wc.species}/{wc.species}.pheno.txt"
    output:
        bed="{bed_output_dir}/{species}/{species}.bed",
        bim="{bed_output_dir}{species}/{species}.bim",
        fam="{bed_output_dir}/{species}/{species}.fam"
    threads: threads
    #conda:
     #   "/home/toriowo/SNAKEMAKE/SexDetermination_snakemake/envs/gemma_env.yaml"
    shell:
        """
        mkdir -p {bed_output_dir}/{wildcards.species}

        out_prefix={bed_output_dir}/{wildcards.species}/{wildcards.species}

        plink --bcf {input.bcf} \
              --double-id \
              --pheno {input.pheno} \
              --make-bed \
              --out {bed_output_dir}/{wildcards.species}/{wildcards.species} \
              --noweb --allow-no-sex
        """

# Rule to run GEMMA association analysis
rule run_gemma:
    input:
        bed=lambda wc: f"{bed_output_dir}/{wc.species}/{wc.species}.bed",
        bim=lambda wc: f"{bed_output_dir}/{wc.species}/{wc.species}.bim",
        fam=lambda wc: f"{bed_output_dir}/{wc.species}/{wc.species}.fam"
    output:
        assoc="{gwas_output_dir}/{species}/output/{species}_gemma.assoc.txt",
        log="{gwas_output_dir}/{species}/output/{species}_gemma.log"
    params:
        in_prefix=lambda wc: f"{bed_output_dir}/{wc.species}/{wc.species}",
        out_prefix=lambda wc: f"{wc.species}_gemma"
    resources:
        mem_mb=8000
    shell:
        """
        mkdir -p {gwas_output_dir}/{wildcards.species}/output

        cd {gwas_output_dir}/{wildcards.species}

        gemma -bfile {params.in_prefix} -lm 2 -o {params.out_prefix} > {gwas_output_dir}/{wildcards.species}/output/{wildcards.species}_gemma.log  2>&1
        """
