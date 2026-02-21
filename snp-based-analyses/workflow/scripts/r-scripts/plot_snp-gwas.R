library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# Define color pairs for each species
# Define color pairs for each species
color_map <- list(
    SeinePhoxinus = c("#9897A9", "#ab9519")
)

process_species_dir <- function(base_dir, species_name, ref_file) {
  species_dir <- file.path(base_dir, species_name)
  output_dir <- file.path(species_dir, "output")
  cat("Processing species:", species_name, "in directory:", output_dir, "\n")

  # Select colors or fallback default
  if (species_name %in% names(color_map)) {
    color1 <- color_map[[species_name]][1]
    color2 <- color_map[[species_name]][2]
  } else {
    color1 <- "#888888"
    color2 <- "#CCCCCC"
    cat("Warning: Species not recognized, using default gray colors.\n")
  }
  cat("Using colors:", color1, color2, "\n")

  gwas_file <- file.path(output_dir, paste0(species_name, "_gemma.assoc.txt"))
  if (!file.exists(gwas_file)) {
    cat("GWAS file not found for species:", species_name, "\n")
    return()
  }

  gwas0 <- read_tsv(gwas_file, show_col_types = FALSE)
  cat("Initial NAs in GWAS data (ps):", sum(is.na(gwas0$ps)), "\n")
  cat("Initial NAs in GWAS data (p_lrt):", sum(is.na(gwas0$p_lrt)), "\n")

  gwas0 <- gwas0 %>%
    rename(scaf = chr) %>%
    mutate(
      ps = ifelse(is.nan(ps), NA, ps),
      p_lrt = ifelse(is.nan(p_lrt), NA, p_lrt)
    )

  cat("NAs after replacing NaN (ps):", sum(is.na(gwas0$ps)), "\n")
  cat("NAs after replacing NaN (p_lrt):", sum(is.na(gwas0$p_lrt)), "\n")

  # Keep only relevant columns & numeric types, drop NAs
  gwas1 <- gwas0 %>%
    select(scaf, ps, n_mis, n_obs, allele0, allele1, af, p_lrt) %>%
    mutate(
      scaf = as.numeric(scaf),
      ps = as.numeric(ps),
      p_lrt = as.numeric(p_lrt)
    ) %>%
    filter(!is.na(ps) & !is.na(p_lrt))

  cat("NAs after filtering (ps):", sum(is.na(gwas1$ps)), "\n")
  cat("NAs after filtering (p_lrt):", sum(is.na(gwas1$p_lrt)), "\n")

  # Top 5% and 1% cutoffs based on p_lrt ascending (smaller p_lrt = stronger)
  cutoff5_n <- max(round(nrow(gwas1) * 0.05), 1)
  cutoff1_n <- max(round(nrow(gwas1) * 0.01), 1)

  gwas_sort_cut5 <- gwas1 %>% arrange(p_lrt) %>% slice_head(n = cutoff5_n)
  perc5 <- max(gwas_sort_cut5$p_lrt)
  cat("Top 5% p-value cutoff:", perc5, "\n")

  gwas_sort_cut1 <- gwas1 %>% arrange(p_lrt) %>% slice_head(n = cutoff1_n)
  perc1 <- max(gwas_sort_cut1$p_lrt)
  cat("Top 1% p-value cutoff:", perc1, "\n")

  # Transform p-values to -log10 scale for plotting
  gwas1 <- gwas1 %>%
    mutate(p_lrt = -log10(p_lrt))

  # Load reference genome info
  ref <- read.table(ref_file, header = FALSE, col.names = c("scaf", "LENGTH", "X1", "X2", "X3")) %>%
    mutate(scaf = as.numeric(gsub('Chr', '', scaf)),
           LENGTH = as.numeric(LENGTH)) %>%
    arrange(scaf) %>%
    mutate(cum_start = cumsum(c(0, LENGTH[-n()])),
           color_group = factor(scaf %% 2))

  # Join GWAS data with reference info
  gwas <- left_join(gwas1, ref, by = "scaf") %>%
    arrange(scaf) %>%
    mutate(Bpos = cum_start + ps)

  cat("NAs after join in Bpos:", sum(is.na(gwas$Bpos)), "\n")

  # Prepare axis breaks (chromosome centers)
  axis_set <- gwas %>%
    group_by(scaf) %>%
    summarize(center = (max(Bpos) + min(Bpos)) / 2, .groups = "drop")

  cols <- c("0" = color1, "1" = color2)

  # Plot GWAS results with top cutoffs indicated as dotted lines
  gwas_plot <- ggplot(gwas, aes(x = Bpos, y = p_lrt, color = color_group)) +
    geom_point(size = 0.5) +
    geom_hline(yintercept = -log10(perc5), linetype = "dotted", size = 1, color = "black") +
    geom_hline(yintercept = -log10(perc1), linetype = "dotted", size = 0.75, color = "black") +
    scale_color_manual(values = cols, guide = "none") +
    scale_x_continuous(labels = axis_set$scaf, breaks = axis_set$center) +
    labs(x = "Chromosome", y = "-log10(p-value)") +
    ggtitle(paste("GWAS outliers for", species_name, "with 5% and 1% cutoff")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  ggsave(file.path(output_dir, paste0(species_name, "_GWAS.tiff")), width = 12, height = 6, plot = gwas_plot)

  # Save data for plotting
  write_tsv(gwas, file = file.path(output_dir, paste0(species_name, "_gwas_plot_values.txt")))

  # Save top 5% GWAS SNPs for combined analysis
  GWAS_out_sort_cut <- gwas0 %>%
    arrange(p_lrt) %>%
    slice_head(n = cutoff5_n) %>%
    select(scaf, ps)

  write_tsv(GWAS_out_sort_cut, file = file.path(output_dir, paste0(species_name, "_GWAS_out_sort_cut.txt")), col_names = FALSE)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if(length(args) < 3) stop("Usage: Rscript script.R <base_dir> <ref_file> <species>")
  
  base_dir <- args[1]
  ref_file <- args[2]
  species_name <- args[3]
  
  process_species_dir(base_dir, species_name, ref_file)
}

main()
