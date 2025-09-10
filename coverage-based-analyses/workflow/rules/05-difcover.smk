# rules/05-difcover.smk

rule difcover:
    input:
        female_bam = lambda wc: os.path.join(config["MERGED_BAM_DIR"], wc.species, f"{wc.species}_Female_merged.bam"),
        male_bam   = lambda wc: os.path.join(config["MERGED_BAM_DIR"], wc.species, f"{wc.species}_Male_merged.bam"),
        cov_tab    = lambda wc: os.path.join(config["OUTPUT_DIR"], "cov_stats", wc.species, f"cov_stats_{wc.species}.tab")
    output:
        ratio_file = os.path.join(config["OUTPUT_DIR"], "difcover", "{species}", "sample1_sample2.ratio_per_w_CC0.txt")
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        mkdir -p $(dirname {output.ratio_file})

        # Parse coverage values
        a=$(awk 'NR==2{{print $2}}' {input.cov_tab})
        A=$(awk 'NR==2{{print $3}}' {input.cov_tab})
        b=$(awk 'NR==2{{print $4}}' {input.cov_tab})
        B=$(awk 'NR==2{{print $5}}' {input.cov_tab})

        # Run pre-step: unionbedcv (optional, skipped here)

        # Run differential coverage
        {config[SCRIPTS_DIR]}/from_unionbed_to_ratio_per_window_CC0 -a $a -A $A -b $b -B $B -v 1000 -l 500 sample1_sample2.unionbedcv > {output.ratio_file}
        """
