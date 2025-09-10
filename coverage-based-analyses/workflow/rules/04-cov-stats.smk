# rules/04-cov-stats.smk

import os

rule cov_stats:
    input:
        female_bam = lambda wc: os.path.join(config["MERGED_BAM_DIR"], wc.species, f"{wc.species}_Female_merged.bam"),
        male_bam   = lambda wc: os.path.join(config["MERGED_BAM_DIR"], wc.species, f"{wc.species}_Male_merged.bam"),
        female_stats = lambda wc: os.path.join(config["MERGED_BAM_DIR"], wc.species, "samtools_stats", f"{wc.species}_Female_merged.stats.txt"),
        male_stats   = lambda wc: os.path.join(config["MERGED_BAM_DIR"], wc.species, "samtools_stats", f"{wc.species}_Male_merged.stats.txt")
    output:
        tab = os.path.join(config["OUTPUT_DIR"], "cov_stats", "{species}", "cov_stats_{species}.tab")
    conda:
        "envs/samtools.yaml"
    run:
        import math
        import subprocess

        os.makedirs(os.path.dirname(output.tab), exist_ok=True)

        def get_min_max(bam):
            depths = subprocess.check_output(f"samtools depth {bam}", shell=True).decode().splitlines()
            covs = [int(line.split()[2]) for line in depths if len(line.split()) == 3]
            return min(covs), max(covs)

        def get_modal_coverage(stats_file):
            modal = 0
            max_freq = 0
            with open(stats_file) as f:
                for line in f:
                    if line.startswith("COV"):
                        fields = line.strip().split("\t")
                        if len(fields) >= 4:
                            cov = int(fields[1])
                            freq = int(fields[3])
                            if freq > max_freq:
                                modal = cov
                                max_freq = freq
            return modal

        # Get all values
        min_f, max_f = get_min_max(input.female_bam)
        min_m, max_m = get_min_max(input.male_bam)
        modal_f = get_modal_coverage(input.female_stats)
        modal_m = get_modal_coverage(input.male_stats)
        cov_ratio = round(modal_m / modal_f, 5) if modal_f != 0 else "NA"

        size_f = os.path.getsize(input.female_bam)
        size_m = os.path.getsize(input.male_bam)
        bam_ratio = round(size_m / size_f, 5) if size_f != 0 else "NA"

        with open(output.tab, "w") as out:
            out.write("SPECIES\tf_min_cov\tf_max_cov\tm_min_cov\tm_max_cov\tcov_ratio\tbam_ratio\n")
            out.write(f"{wildcards.species}\t{min_f}\t{max_f}\t{min_m}\t{max_m}\t{cov_ratio}\t{bam_ratio}\n")
