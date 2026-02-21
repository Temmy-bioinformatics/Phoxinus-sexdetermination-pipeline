04_list_kmers.smk
rule generate_kmers_list:
    input:
        expand(
            f"{config['PROJECT_DIR']}/{config['SPECIES']}/{{sample}}/kmers_with_strand",
            sample=SAMPLES
        )
    output:
        final_list = f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers_list_paths_final.txt",
        dirlist = f"{config['PROJECT_DIR']}/{config['SPECIES']}/dirlist.txt"
    params:
        work_species_dir = f"{config['PROJECT_DIR']}/{config['SPECIES']}"
    shell:
        r"""
        cd {params.work_species_dir}

        # Cleanup intermediate KMC outputs
        find {params.work_species_dir} -type f -name "output_kmc_*" -exec rm -f {{}} +

        # Create dirlist
        ls {params.work_species_dir} | tail -n +2 > {output.dirlist}

        # Build kmers paths
        find {params.work_species_dir} -mindepth 1 -maxdepth 1 -type d | \
        awk -v wd="{params.work_species_dir}" '{{printf "%s/%s/kmers_with_strand\t%s\n", wd, $1, $1}}' \
        > kmers_list_paths_1.txt

        # Remove gz files
        grep -v .gz kmers_list_paths_1.txt > kmers_list_paths_2.txt

        # Keep only dirs in dirlist
        grep -f {output.dirlist} kmers_list_paths_2.txt > kmers_list_paths_clean.txt

        # Fix paths (use | delimiter for sed to avoid issues)
        sed "s|{params.work_species_dir}/||g" kmers_list_paths_clean.txt > {output.final_list}

        # Cleanup temp files
        rm -f kmers_list_paths_clean.txt kmers_list_paths_1.txt kmers_list_paths_2.txt
        """

rule combine_kmers:
    input:
        kmers_path_list = f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers_list_paths_final.txt"
    output:
        f"{config['PROJECT_DIR']}/{config['SPECIES']}/kmers_to_use"
    params:
        work_species_dir = f"{config['PROJECT_DIR']}/{config['SPECIES']}"
    conda:
        "envs/kmer_gwas_envs.yaml"
    shell:
        r"""
        cd {params.work_species_dir}
        {config[KMERS_GWAS_BIN]}/list_kmers_found_in_multiple_samples \
            -l {input.kmers_path_list} \
            -k {config[KMER_LENGTH]} \
            --mac {config[MAC]} \
            -p 0.2 \
            -o {output}
        """
