#06_kmer_gwas.smk
rule kmer_gwas:
    input:
        kmers_table = f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers_table.table",
        kinship     = f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers_table.kinship",
        pheno       = f"{config['PHENO_DIR']}/{config['SPECIES']}.kmer.pheno"
    output:
        assoc     = f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers.assoc",
        assoc_tab = f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers.assoc.tab"
    params:
        work_dir     = f"{config['PROJECT_DIR']}/{config['SPECIES']}",
        table_prefix = f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers_table",
        assoc_prefix = f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers",
        kmers_bin    = config["KMERS_GWAS_BIN"]
    threads: config["THREADS"]
    conda:
        "envs/kmer_gwas_envs.yaml"
    shell:
        r"""
        set -euo pipefail
        cd {params.work_dir}

        {params.kmers_bin}/kmers_table_to_bed \
            -t {params.table_prefix} \
            -k {config[KMER_LENGTH]} \
            -p {input.pheno} \
            --maf {config[MAF]} \
            --mac {config[MAC]} \
            -b {config[BLOCK_SIZE]} \
            -o kmerGWAS_plink

        plink --noweb \
              --bfile kmerGWAS_plink.0 \
              --allow-no-sex \
              --assoc \
              --out {params.assoc_prefix}

        rm -f kmerGWAS_plink.0.*

        awk -v OFS='\t' '{{ $1=$1; print }}' {output.assoc} > {output.assoc_tab}
        """
