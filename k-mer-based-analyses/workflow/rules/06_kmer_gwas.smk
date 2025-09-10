rule gwas:
    input:
        kmers_table = f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers_table.table",
        kinship = f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers_table.kinship",
        pheno = f"{config['PHENO_DIR']}/{config['SPECIES']}.kmer.pheno"
    output:
        assoc = f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers.assoc"
    params:
        work_dir = f"{config['PROJECT_DIR']}/{config['SPECIES']}",
        table_prefix = f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers_table",
        assoc_prefix = f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers"
    threads: config["threads"]
    conda:
        "envs/kmer_gwas_envs.yaml"
    shell:
        """
        cd {params.work_dir}

        {config['kmers_gwas_bin']}/kmers_table_to_bed \
            -t {params.table_prefix} \
            -k {config['kmer_length']} \
            -p {input.pheno} \
            --maf {config['maf']} \
            --mac 1 \
            -b {config['block_size']} \
            -o kmerGWAS_plink

        plink --noweb --bfile kmerGWAS_plink.0 --allow-no-sex --assoc --out {params.assoc_prefix}
        """

rule cleanup_gwas_intermediate:
    input:
        assoc = f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers.assoc"
    output:
        assoc_tab = f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers.assoc.tab"
    params:
        work_dir = f"{config['PROJECT_DIR']}/{config['SPECIES']}"
    threads: config["threads"]
    shell:
        """
        cd {params.work_dir}
        rm -f kmerGWAS_plink.0.*

        awk -v OFS='\\t' '{{ $1=$1; print }}' {input.assoc} > {output.assoc_tab}
        """
