# rules/07-dnacopy.smk
import os

rule dnacopy:
    input:
        log2adj = os.path.join(
            config["OUTPUT_DIR"], "dnacopy", "{species}", "sample1_sample2.ratio_per_w_CC0.log2adj.txt"
        )
    output:
        dnacopy_result = os.path.join(
            config["OUTPUT_DIR"], "dnacopy", "{species}", "sample1_sample2.ratio_per_w_CC0.log2adj.txt.DNAcopyout"
        ),
        dnacopy_plot = os.path.join(
            config["OUTPUT_DIR"], "dnacopy", "{species}", "sample1_sample2.ratio_per_w_CC0.log2adj.txt.pdf"
        )
    params:
        scripts_dir = config["SCRIPTS_DIR"]
    conda:
        "envs/dna-copy-r.yaml"
    shell:
        """
        echo "Running DNAcopy R script..."
        echo "Input file: {input.log2adj}"

        # Run the R script; it will create the two output files in the same directory
        Rscript {params.scripts_dir}/run_DNAcopy_from_bash.R {input.log2adj}

        echo "DNAcopy analysis complete"
        echo "Outputs:"
        echo "- {output.dnacopy_result}"
        echo "- {output.dnacopy_plot}"
        """
