# rules/02-multiqc.smk

rule multiqc:
    input:
        stats_dir=lambda wildcards: os.path.join(
            config["OUTPUT_DIR"], f"samtools_stats/{wildcards.species}"
        )
    output:
        html=os.path.join(config["OUTPUT_DIR"], "multiqc/{species}_multiqc_report.html")
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        mkdir -p $(dirname {output.html})
        multiqc {input.stats_dir} \
            --outdir $(dirname {output.html}) \
            --filename {wildcards.species}_multiqc_report.txt
        """
