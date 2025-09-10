rule add_strand_info:
    input:
        kmc_canon = f"{config['PROJECT_DIR']}/{config['SPECIES']}/{{sample}}/output_kmc_canon.kmc_pre",
        kmc_all = f"{config['PROJECT_DIR']}/{config['SPECIES']}/{{sample}}/output_kmc_all.kmc_pre"
    output:
        f"{config['PROJECT_DIR']}/{config['SPECIES']}/{{sample}}/kmers_with_strand"
    params:
        kmc_canon = f"{config['PROJECT_DIR']}/{config['SPECIES']}/{{sample}}/output_kmc_canon",
        kmc_all = f"{config['PROJECT_DIR']}/{config['SPECIES']}/{{sample}}/output_kmc_all"
    threads: config["THREADS"]
    conda:
        "envs/kmer_gwas_envs.yaml"
    shell:
        """
        cd {config[PROJECT_DIR]}/{config[SPECIES]}/{wildcards.sample}

        {config[TOOLS_DIR]}/kmers_add_strand_information -c {params.kmc_canon} -n {params.kmc_all} -k {config[KMER_LENGTH]} -o {output}
        """
