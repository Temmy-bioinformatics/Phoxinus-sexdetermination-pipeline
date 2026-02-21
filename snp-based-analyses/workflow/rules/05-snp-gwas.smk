import os

project_dir = config["PROJECT_DIR"]
bcf_dir = config["BCF_DIR"]
output_dir = config["OUTPUT_DIR"]
sample_list_dir = config["SAMPLE_DIR"]
threads = config["THREADS"]

bed_output_dir = os.path.join(output_dir, "snp_gwas", "bed_files")
gwas_output_dir = os.path.join(output_dir, "snp_gwas", "gwas")
pheno_output_dir = os.path.join(output_dir, "snp_gwas", "pheno_files")

# ─────────────────────────────────────────────────────────────
# Rule to create phenotype file
# ─────────────────────────────────────────────────────────────
rule create_pheno:
    input:
        male = lambda wc: os.path.join(sample_list_dir, f"{wc.species}_Male.vcflist"),
        female = lambda wc: os.path.join(sample_list_dir, f"{wc.species}_Female.vcflist")
    output:
        pheno = os.path.join(pheno_output_dir, "{species}", "{species}.pheno.txt")
    threads: threads
    shell:
        """
        mkdir -p $(dirname {output.pheno})
        > {output.pheno}
        
        # Process male samples (phenotype = 1)
        while read -r line; do
            # Skip empty lines
            if [ -n "$line" ]; then
                echo -e "${{line}}\t${{line}}\t1" >> {output.pheno}
            fi
        done < {input.male}
        
        # Process female samples (phenotype = 2)
        while read -r line; do
            # Skip empty lines
            if [ -n "$line" ]; then
                echo -e "${{line}}\t${{line}}\t2" >> {output.pheno}
            fi
        done < {input.female}
        """
# ─────────────────────────────────────────────────────────────
# Rule to create PLINK bed/bim/fam files
# ─────────────────────────────────────────────────────────────
rule create_bed:
    input:
        bcf = lambda wc: os.path.join(bcf_dir, wc.species, "Filtered", f"{wc.species}_allChromosomes.filtered.norm.bcf"),
        pheno = lambda wc: os.path.join(pheno_output_dir, wc.species, f"{wc.species}.pheno.txt")
    output:
        bed = os.path.join(bed_output_dir, "{species}", "{species}.bed"),
        bim = os.path.join(bed_output_dir, "{species}", "{species}.bim"),
        fam = os.path.join(bed_output_dir, "{species}", "{species}.fam")
    params:
        out_prefix = lambda wc, output: output.bed[:-4]  # Remove .bed extension
    threads: threads
    conda:
        "envs/gemma_env.yaml"
    shell:
        """
        mkdir -p $(dirname {output.bed})
        
        plink --bcf {input.bcf} \
              --double-id \
              --pheno {input.pheno} \
              --make-bed \
              --out {params.out_prefix} \
              --noweb \
              --allow-no-sex
        """

# ─────────────────────────────────────────────────────────────
# Rule to run GEMMA association analysis
# ─────────────────────────────────────────────────────────────
rule run_gemma:
    input:
        bed = lambda wc: os.path.join(bed_output_dir, wc.species, f"{wc.species}.bed"),
        bim = lambda wc: os.path.join(bed_output_dir, wc.species, f"{wc.species}.bim"),
        fam = lambda wc: os.path.join(bed_output_dir, wc.species, f"{wc.species}.fam")
    output:
        assoc = os.path.join(gwas_output_dir, "{species}", "output", "{species}_gemma.assoc.txt"),
        log = os.path.join(gwas_output_dir, "{species}", "output", "{species}_gemma.log")
    params:
        in_prefix = lambda wc: os.path.join(bed_output_dir, wc.species, wc.species),
        out_prefix = lambda wc: f"{wc.species}_gemma",
        out_dir = lambda wc: os.path.join(gwas_output_dir, wc.species)
    threads: threads
    resources:
        mem_mb = 8000
    conda:
        "envs/gemma_env.yaml"
    shell:
        """
        mkdir -p {params.out_dir}/output
        
        cd {params.out_dir}
        
        gemma -bfile {params.in_prefix} \
              -lm 2 \
              -o {params.out_prefix} > output/{wildcards.species}_gemma.log 2>&1
        """
# ─────────────────────────────────────────────────────────────
# Rule: Plot GWAS results
# ─────────────────────────────────────────────────────────────
rule plot_gwas:
    input:
        assoc=os.path.join(gwas_output_dir, "{species}", "output", "{species}_gemma.assoc.txt"),
        ref=config["ref_file"]
    output:
        plot=os.path.join(gwas_output_dir, "{species}", "output", "{species}_GWAS.tiff"),
        plot_values=os.path.join(gwas_output_dir, "{species}", "output", "{species}_gwas_plot_values.txt"),
        top_snps=os.path.join(gwas_output_dir, "{species}", "output", "{species}_GWAS_out_sort_cut.txt")
    params:
        script=config["gwas_rscript"],
        gwas_base_dir=gwas_output_dir
    conda:
        "envs/r-plots_env.yaml"
    shell:
        """
        Rscript {params.script} {params.gwas_base_dir} {input.ref} {wildcards.species}
        """

# ─────────────────────────────────────────────────────────────
# Rule: Sort top GWAS SNPs by coordinates
# ─────────────────────────────────────────────────────────────
rule sort_top_gwas_snps:
    input:
        snps=os.path.join(gwas_output_dir, "{species}", "output", "{species}_GWAS_out_sort_cut.txt")
    output:
        sorted=os.path.join(gwas_output_dir, "{species}", "output", "{species}_coordinated_sorted_GWAS_cut.txt")
    shell:
        """
        sort -k1,1 -k2n {input.snps} > {output.sorted}
        """

# ─────────────────────────────────────────────────────────────
# Rule: Count top GWAS SNPs per 10kb window
# ─────────────────────────────────────────────────────────────
rule count_gwas_windows:
    input:
        sorted=os.path.join(gwas_output_dir, "{species}", "output", "{species}_coordinated_sorted_GWAS_cut.txt")
    output:
        counts=os.path.join(gwas_output_dir, "{species}", "output", "{species}_GWAS_window_count.txt")
    params:
        script=config["count_window"],
        window_size=10000
    shell:
        """
        python {params.script} {input.sorted} {output.counts} {params.window_size}
        """