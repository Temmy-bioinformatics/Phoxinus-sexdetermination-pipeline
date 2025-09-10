# rules/06-dnacopy.smk

rule dnacopy_input:
    input:
        ratio_file = os.path.join(config["OUTPUT_DIR"], "difcover", "{species}", "sample1_sample2.ratio_per_w_CC0.txt"),
        cov_tab    = os.path.join(config["OUTPUT_DIR"], "cov_stats", "{species}", "cov_stats_{species}.tab")
    output:
        log2adj = os.path.join(config["OUTPUT_DIR"], "dnacopy", "{species}", "sample1_sample2.ratio_per_w_CC0.log2adj.txt")
    conda:
        "envs/deeptools.yaml"  # Just to get Python and awk; feel free to replace
    shell:
        """
        mkdir -p $(dirname {output.log2adj})

        AC=$(awk 'NR==2{{print $7}}' {input.cov_tab})

        echo -e "scaffold\twindow_start\tAC*log2(ratio)" > {output.log2adj}
        awk -v AC="$AC" 'BEGIN{{OFS="\\t"}} $2 != "window_start" {{print $1, $2, log(AC * $7)/log(2)}}' {input.ratio_file} >> {output.log2adj}
        """
