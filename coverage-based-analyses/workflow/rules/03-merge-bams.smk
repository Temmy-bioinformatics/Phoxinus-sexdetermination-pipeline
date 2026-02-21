# rules/03-merge-bams.smk

rule merge_bams:
    input:
        bam_list=lambda wildcards: os.path.join(
            config["PROJECT_DIR"],
            "resources/bamlists",
            f"{wildcards.species}_{wildcards.sex}_Bamfiles_Metadata.txt"
        )
    output:
        merged=os.path.join(config["MERGED_BAM_DIR"], "{species}", "{species}_{sex}_merged.bam"),
        index=os.path.join(config["MERGED_BAM_DIR"], "{species}", "{species}_{sex}_merged.bam.bai"),
        stats=os.path.join(config["MERGED_BAM_DIR"], "{species}", "samtools_stats", "{species}_{sex}_merged.stats.txt")
    conda:
        "envs/samtools.yaml"
    threads: 8
    shell:
        """
        mkdir -p $(dirname {output.merged})
        mkdir -p $(dirname {output.stats})

        samtools merge -@ {threads} -b {input.bam_list} {output.merged}
        samtools index {output.merged}
        samtools stats -@ {threads} {output.merged} > {output.stats}
        """
