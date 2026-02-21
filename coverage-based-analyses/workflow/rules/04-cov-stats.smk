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
        import subprocess
        
        os.makedirs(os.path.dirname(output.tab), exist_ok=True)
        
        def get_coverage_stats(bam):
            """Use samtools coverage for much faster min/max coverage calculation"""
            result = subprocess.check_output(
                f"samtools coverage {bam}", 
                shell=True
            ).decode().splitlines()
            
            # Skip header, process data lines
            coverages = []
            for line in result[1:]:  # Skip header
                fields = line.split('\t')
                if len(fields) >= 7:
                    meandepth = float(fields[6])  # meandepth column
                    coverages.append(meandepth)
            
            if coverages:
                return min(coverages), max(coverages), sum(coverages) / len(coverages)
            return 0, 0, 0
        
        def get_modal_coverage(stats_file):
            """Extract modal coverage from samtools stats COV section"""
            modal = 0
            max_freq = 0
            with open(stats_file) as f:
                for line in f:
                    if line.startswith("COV"):
                        fields = line.strip().split("\t")
                        if len(fields) >= 4:
                            try:
                                cov = int(fields[2])  # Coverage depth
                                freq = int(fields[3])  # Frequency
                                if freq > max_freq:
                                    modal = cov
                                    max_freq = freq
                            except ValueError:
                                continue
            return modal
        
        # Get coverage stats
        min_f, max_f, mean_f = get_coverage_stats(input.female_bam)
        min_m, max_m, mean_m = get_coverage_stats(input.male_bam)
        
        # Get modal coverage from stats files
        modal_f = get_modal_coverage(input.female_stats)
        modal_m = get_modal_coverage(input.male_stats)
        
        cov_ratio = round(modal_m / modal_f, 5) if modal_f != 0 else "NA"
        
        # BAM file size ratio
        size_f = os.path.getsize(input.female_bam)
        size_m = os.path.getsize(input.male_bam)
        bam_ratio = round(size_m / size_f, 5) if size_f != 0 else "NA"
        
        with open(output.tab, "w") as out:
            out.write("SPECIES\tf_min_cov\tf_max_cov\tm_min_cov\tm_max_cov\tf_modal\tm_modal\tcov_ratio\tbam_ratio\n")
            out.write(f"{wildcards.species}\t{min_f:.2f}\t{max_f:.2f}\t{min_m:.2f}\t{max_m:.2f}\t{modal_f}\t{modal_m}\t{cov_ratio}\t{bam_ratio}\n")