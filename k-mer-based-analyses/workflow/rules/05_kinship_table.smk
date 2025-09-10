rule build_kmers_table:
    input:
        kmers_to_use = f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers_to_use",
        kmers_path_list = f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers_list_paths_final.txt",
        kmers_with_strand = expand(
            f"{config['PROJECT_DIR']}/{config['SPECIES']}/{{sample}}/kmers_with_strand",
            sample=SAMPLES
        )
    output:
        f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers_table.table"
    params:
        work_species_dir = f"{config['PROJECT_DIR']}/{config['SPECIES']}",
        table_prefix = f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers_table"
    conda:
        "envs/kmer_gwas_envs.yaml"
    shell:
        """
        cd {params.work_species_dir}
        {config[KMERS_GWAS_BIN]}/build_kmers_table \
            -l {input.kmers_path_list} \
            -k {config[KMER_LENGTH]} \
            -a {input.kmers_to_use} \
            -o {params.table_prefix}
        """

rule kinship:
    input:
        f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers_table.table"
    output:
        f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers_table.kinship"
    params:
        work_species_dir = f"{config['PROJECT_DIR']}/{config['SPECIES']}",
        table_prefix = f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers_table"
    conda:
        "envs/kmer_gwas_envs.yaml"
    shell:
        """
        cd {params.work_species_dir}
        {config[KMERS_GWAS_BIN]}/emma_kinship_kmers \
            -t {params.table_prefix} \
            -k {config[KMER_LENGTH]} \
            --maf {config[MAF]} > {output}
        """
