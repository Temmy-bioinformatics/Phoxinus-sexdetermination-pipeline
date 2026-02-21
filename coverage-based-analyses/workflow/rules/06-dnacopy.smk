# rules/06-dnacopy.smk
rule dnacopy_input:
    input:
        ratio_dir = os.path.join(config["OUTPUT_DIR"], "difcover", "{species}"),
        cov_tab = os.path.join(config["OUTPUT_DIR"], "cov_stats", "{species}", "cov_stats_{species}.tab")
    output:
        log2adj = os.path.join(config["OUTPUT_DIR"], "dnacopy", "{species}", "sample1_sample2.ratio_per_w_CC0.log2adj.txt")
    params:
        ratio_pattern = lambda wc: os.path.join(config["OUTPUT_DIR"], "difcover", wc.species, "sample1_sample2.ratio_per_w_CC0_a*")
    conda:
        "envs/deeptools.yaml"
    shell:
        """
        mkdir -p $(dirname {output.log2adj})
        
        # Find the ratio file (it has parameters in the filename)
        ratio_file=$(ls {params.ratio_pattern} 2>/dev/null | head -1)
        
        if [ -z "$ratio_file" ] || [ ! -f "$ratio_file" ]; then
            echo "Error: Could not find ratio file matching {params.ratio_pattern}"
            ls -lh {input.ratio_dir}/
            exit 1
        fi
        
        echo "Using ratio file: $ratio_file"
        
        # Extract AC value from coverage stats
        AC=$(awk 'NR==2{{print $7}}' {input.cov_tab})
        echo "AC ratio: $AC"
        
        # Create header
        echo -e "scaffold\\twindow_start\\t${{AC}}*log2(ratio)" > {output.log2adj}
        
        # Process the ratio file - skip parameter lines and header
        awk -v AC="$AC" '
        BEGIN {{
            OFS="\\t"
        }}
        /^Parameters:/ {{ next }}
        /^[a-z]=/ {{ next }}
        /^scaffold/ {{ next }}
        NF >= 7 && $7 ~ /^[0-9]+\\.?[0-9]*$/ && $7 > 0 {{
            result = log(AC * $7) / log(2)
            print $1, $2, result
        }}
        ' "$ratio_file" >> {output.log2adj}
        
        # Verify output
        line_count=$(wc -l < {output.log2adj})
        echo "Log2adj file created with $line_count lines (including header)"
        
        if [ "$line_count" -lt 2 ]; then
            echo "WARNING: No data lines found"
            exit 1
        fi
        
        echo "Successfully processed ratio file"
        """