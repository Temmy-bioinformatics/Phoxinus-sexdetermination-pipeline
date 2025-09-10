rule copy_male_fastq:
    input:
        fastq1 = lambda wildcards: f"{config['SOURCE_DIRECTORY']}/{config['SPECIES']}/Male/{wildcards.sample}_1.fq.gz",
        fastq2 = lambda wildcards: f"{config['SOURCE_DIRECTORY']}/{config['SPECIES']}/Male/{wildcards.sample}_2.fq.gz"
    output:
        fastq1 = f"{config['PROJECT_DIR']}/{config['SPECIES']}/{{sample}}/{{sample}}_1.fq.gz",
        fastq2 = f"{config['PROJECT_DIR']}/{config['SPECIES']}/{{sample}}/{{sample}}_2.fq.gz"
    threads: 4
    run:
        if os.path.exists(input.fastq1) and os.path.exists(input.fastq2):
            shell("cp {input.fastq1} {output.fastq1}")
            shell("cp {input.fastq2} {output.fastq2}")
        else:
            print(f"ERROR: Required male files for {wildcards.sample} not found.")
            raise FileNotFoundError(f"Missing male FASTQ files for {wildcards.sample}.")

rule copy_female_fastq:
    input:
        fastq1 = lambda wildcards: f"{config['SOURCE_DIRECTORY']}/{config['SPECIES']}/Female/{wildcards.sample}_1.fq.gz",
        fastq2 = lambda wildcards: f"{config['SOURCE_DIRECTORY']}/{config['SPECIES']}/Female/{wildcards.sample}_2.fq.gz"
    output:
        fastq1 = f"{config['PROJECT_DIR']}/{config['SPECIES']}/{{sample}}/{{sample}}_1.fq.gz",
        fastq2 = f"{config['PROJECT_DIR']}/{config['SPECIES']}/{{sample}}/{{sample}}_2.fq.gz"
    threads: 4
    run:
        if os.path.exists(input.fastq1) and os.path.exists(input.fastq2):
            shell("cp {input.fastq1} {output.fastq1}")
            shell("cp {input.fastq2} {output.fastq2}")
        else:
            print(f"ERROR: Required female files for {wildcards.sample} not found.")
            raise FileNotFoundError(f"Missing female FASTQ files for {wildcards.sample}.")
