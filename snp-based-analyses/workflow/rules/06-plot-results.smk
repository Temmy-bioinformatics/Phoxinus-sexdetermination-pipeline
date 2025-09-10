rule plot_fst:
    input:
        fst = "results/fst/{species}/{species}_Fst_10k_window.txt",
        ref = "resources/reference_genome.fai"
    output:
        plot = "results/plotting/{species}/{species}_Fst_10k_plot.tiff",
        thresholds = "results/plotting/{species}/{species}_Fst_10k_thresholds.txt"
    params:
        script = lambda wildcards: os.path.join(config["RSCRIPT_DIR"], "plot_fst.r")
    conda:
        "workflow/envs/r_environment.yaml"
    shell:
        """
        Rscript {params.script} {input.fst} {input.ref} {wildcards.species} {output.plot} {output.thresholds}
        """

rule plot_nucleotide_diversity:
    input:
        nucdiv = "results/nucleotide_diversity/{species}/{species}_NucDiv_10kb_sex_differences.txt"
    output:
        plot = "results/plotting/{species}/{species}_NucDiv_10kb_sex_differences_plot.tiff"
    params:
        script = lambda wildcards: os.path.join(config["RSCRIPT_DIR"], "plot_nucleotide_diversity.r")
    conda:
        "workflow/envs/r_environment.yaml"
    shell:
        """
        Rscript {params.script} {input.nucdiv} {wildcards.species} {output.plot}
        """

rule plot_snpdensity:
    input:
        snpdensity = "results/snpdensity/{species}/plot/{species}_SNPdensity_SexFindR.txt",
        ref = "resources/reference_genome.fai"
    output:
        plot = "results/plotting/{species}/{species}_SNPdensity_plot.pdf",
        significant = "results/plotting/{species}/{species}_significant_SNPdensity_SexFindR.txt"
    params:
        script = lambda wildcards: os.path.join(config["RSCRIPT_DIR"], "plot_snpdensity.r")
    conda:
        "workflow/envs/r_environment.yaml"
    shell:
        """
        Rscript {params.script} {input.snpdensity} {input.ref} {wildcards.species} {output.plot} {output.significant}
        """

rule plot_gwas:
    input:
        gwas = "results/gwas/{species}/output/{species}_gwas_plot_values.txt",
        ref = "resources/reference_genome.fai"
    output:
        plot = "results/plotting/{species}/{species}_GWAS.tiff",
        top_snps = "results/plotting/{species}/{species}_GWAS_out_sort_cut.txt"
    params:
        script = lambda wildcards: os.path.join(config["RSCRIPT_DIR"], "plot_gwas.r")
    conda:
        "workflow/envs/r_environment.yaml"
    shell:
        """
        Rscript {params.script} {input.gwas} {input.ref} {wildcards.species} {output.plot} {output.top_snps}
        """
