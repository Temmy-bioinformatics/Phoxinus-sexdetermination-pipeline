#08_blast_kmers.smk
rule process_top_kmers:
    input:
        most_significant_assoc = f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers.assoc.tab",
        most_assoc_kmers = f"{config['PROJECT_DIR']}/{config['SPECIES']}/most_assoc_kmers.list"
    output:
        assoc_kmers_presab = f"{config['PROJECT_DIR']}/{config['SPECIES']}/assoc_kmers_presab.txt",
        plink_abyss_input = f"{config['PROJECT_DIR']}/{config['SPECIES']}/plink_abyss_input.txt",
        abyss_output = f"{config['PROJECT_DIR']}/{config['SPECIES']}/plink_abyss_output.txt",
        blast_output = f"{config['PROJECT_DIR']}/{config['SPECIES']}/blast/kmers_blast.out",
        mean_contig_len = f"{config['PROJECT_DIR']}/{config['SPECIES']}/blast/mean_contig_len"
    params:
        work_dir = f"{config['PROJECT_DIR']}/{config['SPECIES']}",
        table_prefix = f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers_table",
        reference = config['REFERENCE'],
        kmers_gwas_bin = config['KMERS_GWAS_BIN']
    threads: config["THREADS"]
    conda:
        "envs/kmer_gwas_envs.yaml"
    shell:
        """
        cd {params.work_dir}

        {params.kmers_gwas_bin}/filter_kmers -t {params.table_prefix} -k {input.most_assoc_kmers} -o {output.assoc_kmers_presab}

        python /home/toriowo/PROJECTS/sex-determination-pipeline/code/job-scripts/k-mer-based-analyses/workflow/scripts/plink_to_abyss_kmers.py \
            {input.most_significant_assoc} {output.plink_abyss_input}

        ABYSS -k25 -c0 -e0 {output.plink_abyss_input} -o {output.abyss_output}

        makeblastdb -in {params.reference} -dbtype nucl

        blastn -query {output.abyss_output} -db {params.reference} -outfmt 6 -out {output.blast_output}

        cat {output.abyss_output} | grep ">*" | cut -f 2 -d " " -s > {output.mean_contig_len}
        """
