import os

rule samtools_stats:
    input:
        bam_metadata=lambda wildcards: os.path.join(
            config["PROJECT_DIR"],
            "resources/bamlists",
            f"{wildcards.species}_{wildcards.sex}_Bamfiles_Metadata.txt"
        )
    output:
        stats_dir=directory(os.path.join(config["OUTPUT_DIR"], "samtools_stats/{species}/{sex}"))
    conda:
        "envs/samtools.yaml"
    run:
        species = wildcards.species
        sex = wildcards.sex
        metadata = input.bam_metadata
        outdir = output.stats_dir

        os.makedirs(outdir, exist_ok=True)
        logdir = os.path.join(config["OUTPUT_DIR"], "logs/samtools_stats")
        os.makedirs(logdir, exist_ok=True)

        with open(metadata, "r") as f:
            for line in f:
                bam = line.strip()
                bamfile = os.path.basename(bam)
                out_file = os.path.join(outdir, f"{bamfile}.stats.txt")
                err_file = os.path.join(logdir, f"{bamfile}.err")

                shell(f"""
                    echo "[{species} - {sex}] Running samtools stats for {bamfile}"
                    samtools stats {bam} > {out_file} 2> {err_file}
                """)
