# 01_files_prep.smk

rule copy_fastq:
    input:
        fastq1 = lambda wildcards: f"{config['SOURCE_DIRECTORY']}/{species_list[0]}/{get_sample_gender(wildcards.sample)}/{wildcards.sample}_1.fq.gz",
        fastq2 = lambda wildcards: f"{config['SOURCE_DIRECTORY']}/{species_list[0]}/{get_sample_gender(wildcards.sample)}/{wildcards.sample}_2.fq.gz"
    output:
        fastq1 = f"{config['PROJECT_DIR']}/{species_list[0]}/{{sample}}/{{sample}}_1.fq.gz",
        fastq2 = f"{config['PROJECT_DIR']}/{species_list[0]}/{{sample}}/{{sample}}_2.fq.gz"
    threads: 4
    shell:
        """
        cp {input.fastq1} {output.fastq1}
        cp {input.fastq2} {output.fastq2}
        """