# Load required libraries
library(tidyverse)
library(patchwork)
library(ggthemes)

# Main function
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Check for correct number of arguments
  if (length(args) != 11) {
    stop(
      "Usage: Rscript sexfindr-combined-plot.R ",
      "<snp_density_file> <gwas_window_file> <fst_window_file> <pi_diff_file> ",
      "<ref_file> <species_name> ",
      "<output_candidates_all> ",
      "<output_candidates_gwas_pi_fst> ",
      "<output_candidates_snp_fst_pi> ",
      "<output_correlations> ",
      "<output_plot>"
    )
  }

  # Parse arguments
  snp_file <- args[1]
  gwas_file <- args[2]
  fst_file <- args[3]
  pi_file <- args[4]
  ref_file <- args[5]
  species <- args[6]
  output_candidates_all <- args[7]
  output_candidates_gwas_pi_fst <- args[8]
  output_candidates_snp_fst_pi <- args[9]
  output_correlations <- args[10]
  output_plot <- args[11]
  
  # Log processing information
  cat("================================================\n")
  cat("Integrated Analysis for:", species, "\n")
  cat("================================================\n")
  cat("SNP density file:     ", snp_file, "\n")
  cat("GWAS window file:     ", gwas_file, "\n")
  cat("FST window file:      ", fst_file, "\n")
  cat("Pi differences file:  ", pi_file, "\n")
  cat("Reference file:       ", ref_file, "\n")
  cat("Output candidates:    ", output_candidates_all, "\n")
  cat("Output correlations:  ", output_correlations, "\n")
  cat("Output plot:          ", output_plot, "\n")
  cat("================================================\n\n")
  
  # Check if input files exist
  if (!file.exists(snp_file)) stop("SNP density file not found: ", snp_file)
  if (!file.exists(gwas_file)) stop("GWAS window file not found: ", gwas_file)
  if (!file.exists(fst_file)) stop("FST window file not found: ", fst_file)
  if (!file.exists(pi_file)) stop("Pi file not found: ", pi_file)
  if (!file.exists(ref_file)) stop("Reference file not found: ", ref_file)
  
  # ------------------------
  # LOAD AND PROCESS DATA
  # ------------------------
  
  cat("Loading SNP density data...\n")
  SNP <- read_delim(snp_file, show_col_types = FALSE) %>%
    select(scaf, base, mean_MvF_dif, Pvalue) %>%
    mutate(
      scaf = str_remove(scaf, "Chr"),
      base = base + 10000
    ) %>%
    filter(Pvalue <= 0.05)
  
  # Rank SNP density windows
  SNP_filter_arrange <- SNP %>%
    arrange(Pvalue, desc(abs(mean_MvF_dif))) %>%
    mutate(SNPdensity_rank = row_number(),
           scaf = as.numeric(scaf))
  
  cat("  - Loaded", nrow(SNP), "significant SNP density windows\n")
  
  cat("Loading GWAS window counts...\n")
  GWAS_window_count <- read_delim(gwas_file, show_col_types = FALSE) %>%
    arrange(desc(count)) %>%
    mutate(GWAScount_rank = row_number(),
           scaf = as.numeric(str_remove(scaf, "Chr")))
  
  cat("  - Loaded", nrow(GWAS_window_count), "GWAS windows\n")
  
  cat("Loading FST window counts...\n")
  Fst_window_count <- read_delim(fst_file, show_col_types = FALSE) %>%
    arrange(desc(count)) %>%
    mutate(Fstcount_rank = row_number(),
           scaf = as.numeric(str_remove(scaf, "Chr")))
  
  cat("  - Loaded", nrow(Fst_window_count), "FST windows\n")
  
  cat("Loading nucleotide diversity data...\n")
  pi <- read_delim(pi_file, show_col_types = FALSE) %>%
    select(scaf, BIN_START, pi_dif_M_minus_F) %>%
    arrange(desc(abs(pi_dif_M_minus_F))) %>%
    mutate(pi_rank = row_number(),
           base = BIN_START + 9999,
           scaf = as.numeric(scaf))
  
  cat("  - Loaded", nrow(pi), "pi windows\n")
  
  # Load reference chromosome lengths and calculate cumulative starts
  cat("Loading reference genome...\n")
  ref <- read.table(ref_file, header = FALSE, col.names = c("scaf", "LENGTH", "X1", "X2", "X3")) %>%
    mutate(scaf = as.numeric(str_remove(scaf, "Chr"))) %>%
    arrange(scaf) %>%
    mutate(cum_start = cumsum(c(0, LENGTH[-length(LENGTH)])))
  
  cat("  - Found", nrow(ref), "chromosomes\n")
  
  # Define color map for different species
  color_map <- list(
    SeinePhoxinus = c("#9897A9", "#ab9519"),
    RhinePhoxinus = c("#9897A9", "#66A61E")
  )
  
  # Use species-specific colors or default
  color_pairs <- if (species %in% names(color_map)) {
    color_map[[species]]
  } else {
    c("#888888", "#CCCCCC")
  }
  
  cat("Using colors:", color_pairs[1], ",", color_pairs[2], "\n\n")
  
  # Assign color groups alternating by chromosome
  ref <- ref %>%
    mutate(color_group = color_pairs[(scaf %% 2) + 1])
  
  # Midpoints for x-axis labeling
  axis_set <- ref %>%
    group_by(scaf) %>%
    summarize(center = cum_start + LENGTH / 2, .groups = "drop")
  
  # ------------------------
  # JOIN DATASETS
  # ------------------------
  
  cat("Joining all datasets...\n")
  Fst_window_count_chr <- left_join(Fst_window_count, ref, by = "scaf")
  GWAS_window_count_chr <- left_join(GWAS_window_count, ref, by = "scaf")
  SNP_filter_arrange_chr <- left_join(SNP_filter_arrange, ref, by = "scaf")
  
  SNP_Fst <- full_join(SNP_filter_arrange_chr, Fst_window_count_chr, by = c("scaf", "base"))
  SNP_Fst_GWAS <- left_join(SNP_Fst, GWAS_window_count, by = c("scaf", "base"))
  
  all <- full_join(SNP_Fst_GWAS, pi, by = c("scaf", "base")) %>%
    select(-count.x, -count.y) %>%
    mutate(total_rank = GWAScount_rank + SNPdensity_rank + Fstcount_rank + pi_rank)
  
  # ------------------------
  # SPEARMAN CORRELATION OF RANKINGS
  # ------------------------
  
  cat("Calculating Spearman correlations...\n")
  correlations <- c(
    "SNP_density vs GWAS_count" = cor(all$SNPdensity_rank, all$GWAScount_rank, method = "spearman", use = "complete.obs"),
    "SNP_density vs Fst_count" = cor(all$SNPdensity_rank, all$Fstcount_rank, method = "spearman", use = "complete.obs"),
    "SNP_density vs pi" = cor(all$SNPdensity_rank, all$pi_rank, method = "spearman", use = "complete.obs"),
    "GWAS_count vs Fst_count" = cor(all$GWAScount_rank, all$Fstcount_rank, method = "spearman", use = "complete.obs"),
    "GWAS_count vs pi" = cor(all$GWAScount_rank, all$pi_rank, method = "spearman", use = "complete.obs"),
    "Fst_count vs pi" = cor(all$Fstcount_rank, all$pi_rank, method = "spearman", use = "complete.obs")
  )
  
  # Write sorted correlations to file
  cat("Writing correlation results...\n")
  sink(output_correlations)
  cat("Ranked Spearman's Correlations (from strongest to weakest):\n\n")
  sorted_correlations <- sort(abs(correlations), decreasing = TRUE)
  for (name in names(sorted_correlations)) {
    cat(name, " - rho = ", correlations[name], "\n")
  }
  sink()
  cat("  ✓ Correlations saved to:", output_correlations, "\n")
  
  # ------------------------
  # FILTER TOP WINDOWS BY RANK (top 100)
  # ------------------------
  
  cat("Filtering top-ranked windows...\n")
  top_all <- all %>%
    filter(GWAScount_rank <= 100,
           Fstcount_rank <= 100,
           pi_rank <= 100,
           SNPdensity_rank <= 100)
  
  cat("  - Found", nrow(top_all), "windows in top 100 across all methods\n")
  
  # 3-way filters
  top_rank_SFP <- all %>% 
    filter(SNPdensity_rank <= 100, Fstcount_rank <= 100, pi_rank <= 100)
  
  top_rank_GFP <- all %>% 
    filter(GWAScount_rank <= 100, Fstcount_rank <= 100, pi_rank <= 100)
  
  cat("  - Found", nrow(top_rank_GFP), "windows in top 100 for GWAS+FST+Pi\n")
  cat("  - Found", nrow(top_rank_SFP), "windows in top 100 for SNP+FST+Pi\n")
  
  # Save filtered datasets
  cat("Writing candidate window files...\n")
  write_tsv(top_all %>% select(-total_rank), output_candidates_all)
  cat("  ✓ All methods:", output_candidates_all, "\n")
  
  write_tsv(top_rank_GFP %>% select(-total_rank), output_candidates_gwas_pi_fst)
  cat("  ✓ GWAS+Pi+FST:", output_candidates_gwas_pi_fst, "\n")
  
  write_tsv(top_rank_SFP %>% select(-total_rank), output_candidates_snp_fst_pi)
  cat("  ✓ SNP+FST+Pi:", output_candidates_snp_fst_pi, "\n")
  
  # ------------------------
  # PREPARE WINDOW POSITIONS FOR PLOTTING
  # ------------------------
  
  add_Wpos <- function(df) {
    df %>%
      mutate(Wpos = cum_start + base)
  }
  
  Fst_window_count_chr <- add_Wpos(Fst_window_count_chr)
  GWAS_window_count_chr <- add_Wpos(GWAS_window_count_chr)
  SNP_filter_arrange_chr <- add_Wpos(SNP_filter_arrange_chr)
  pi_chr <- add_Wpos(left_join(pi, ref, by = "scaf"))
  
  # FIX: Join top_all with ref to get cum_start column before adding Wpos
  top_all_with_ref <- top_all %>%
    left_join(select(ref, scaf, cum_start), by = "scaf")
  
  top_points <- add_Wpos(top_all_with_ref)
  top_window <- top_points$Wpos
  
  # ------------------------
  # PLOTTING
  # ------------------------
  
  cat("Generating integrated plots...\n")
  
  # FST window counts
  plot_Fst <- Fst_window_count_chr %>%
    ggplot(aes(x = Wpos, y = count, color = color_group)) +
    geom_point(size = 2) +
    geom_vline(xintercept = top_window, linetype = "dotted", color = "red", linewidth = 0.5, alpha = 0.5) +
    scale_x_continuous(breaks = axis_set$center, labels = axis_set$scaf) +
    scale_color_manual(values = color_pairs) +
    labs(
      subtitle = expression("A. Intersex F"[ST]),
      y = expression("F"[ST] * " outliers per 10k window"),
      x = NULL
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.subtitle = element_text(face = "bold")
    )
  
  # GWAS window counts
  plot_GWAS <- GWAS_window_count_chr %>%
    ggplot(aes(x = Wpos, y = count, color = color_group)) +
    geom_point(size = 2) +
    geom_vline(xintercept = top_window, linetype = "dotted", color = "red", linewidth = 0.5, alpha = 0.5) +
    scale_x_continuous(breaks = axis_set$center, labels = axis_set$scaf) +
    scale_color_manual(values = color_pairs) +
    labs(
      subtitle = "B. GWAS",
      y = "GWAS outliers per 10k window",
      x = NULL
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.subtitle = element_text(face = "bold")
    )
  
  # SNP density
  plot_SNP <- SNP_filter_arrange_chr %>%
    ggplot(aes(x = Wpos, y = mean_MvF_dif, color = color_group)) +
    geom_point(size = 2) +
    geom_vline(xintercept = top_window, linetype = "dotted", color = "red", linewidth = 0.5, alpha = 0.5) +
    scale_x_continuous(breaks = axis_set$center, labels = axis_set$scaf) +
    scale_color_manual(values = color_pairs) +
    labs(
      subtitle = "C. SNP Density",
      y = "Mean difference (M - F)",
      x = NULL
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.subtitle = element_text(face = "bold")
    )
  
  # Nucleotide diversity
  plot_Pi <- pi_chr %>%
    ggplot(aes(x = Wpos, y = pi_dif_M_minus_F, color = color_group)) +
    geom_point(size = 2) +
    geom_vline(xintercept = top_window, linetype = "dotted", color = "red", linewidth = 0.5, alpha = 0.5) +
    scale_x_continuous(breaks = axis_set$center, labels = axis_set$scaf) +
    scale_color_manual(values = color_pairs) +
    labs(
      subtitle = expression("D. Nucleotide Diversity (" * pi * ")"),
      y = expression(pi * " difference (M - F)"),
      x = "Chromosome"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      plot.subtitle = element_text(face = "bold")
    )
  
  # Combine plots
  combined_plot <- plot_Fst / plot_GWAS / plot_SNP / plot_Pi +
    plot_annotation(
      title = paste("Integrated Analysis for", species),
      subtitle = "Red dotted lines indicate windows ranking in top 100 across all four methods",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    )
  
  # Save plot
  ggsave(
    filename = output_plot,
    plot = combined_plot,
    width = 14,
    height = 12,
    dpi = 300
  )
  
  cat("  ✓ Plot saved to:", output_plot, "\n")
  
  cat("\n================================================\n")
  cat("Integrated analysis completed successfully!\n")
  cat("================================================\n")
}

# Run main
main()