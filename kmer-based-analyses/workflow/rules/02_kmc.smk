#02_kmc.smk
rule kmc:
    input:
        fastq1 = f"{config['PROJECT_DIR']}/{config['SPECIES']}/{{sample}}/{{sample}}_1.fq.gz",
        fastq2 = f"{config['PROJECT_DIR']}/{config['SPECIES']}/{{sample}}/{{sample}}_2.fq.gz"
    output:
        kmc_canon_pre = f"{config['PROJECT_DIR']}/{config['SPECIES']}/{{sample}}/output_kmc_canon.kmc_pre",
        kmc_canon_suf = f"{config['PROJECT_DIR']}/{config['SPECIES']}/{{sample}}/output_kmc_canon.kmc_suf",
        kmc_all_pre   = f"{config['PROJECT_DIR']}/{config['SPECIES']}/{{sample}}/output_kmc_all.kmc_pre",
        kmc_all_suf   = f"{config['PROJECT_DIR']}/{config['SPECIES']}/{{sample}}/output_kmc_all.kmc_suf"
    threads: config["THREADS"]
    log:
        f"{config['PROJECT_DIR']}/{config['SPECIES']}/{{sample}}/kmc.log"
    shell:
        r"""
        set -euo pipefail

        cd "{config[PROJECT_DIR]}/{config[SPECIES]}/{wildcards.sample}"

        # Create a file listing input FASTQs
        printf "%s\n%s\n" "{input.fastq1}" "{input.fastq2}" > input_files.txt

        # Run KMC for canonical kmers
        kmc -t{threads} -k{config[KMER_LENGTH]} -ci2 @input_files.txt output_kmc_canon ./ 1> kmc_canon.1 2> kmc_canon.2

        # Run KMC for all kmers (including non-canonical)
        kmc -t{threads} -k{config[KMER_LENGTH]} -ci0 -b @input_files.txt output_kmc_all ./ 1> kmc_all.1 2> kmc_all.2
        """
