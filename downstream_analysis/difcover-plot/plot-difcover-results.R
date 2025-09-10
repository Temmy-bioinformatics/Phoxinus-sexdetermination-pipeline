# Load required libraries
library(tidyverse)
library(ggpubr)
library(ggforce)

# Define color pairs for each species
color_map <- list(
  "WeserMorella"   = c("#66A61E", "darkgrey"),
  "RhineCsikii"    = c("#CCBFFF", "#654CFF"),
  "RhinePhoxinus"  = c("#A5EDFF", "#5B84B1FF"),
  "DanubeCsikii"   = c("#D95E02", "darkgrey")
)

# Reference file (scaffold lengths)
ref_file <- "Final_Phoxy_1_rearranged.fasta.fai"

# Function to plot DifCover results for a species and haplotype
plot_difcover_results <- function(species_name, haplotype) {
  message("Processing ", species_name, " - ", haplotype)

  working_dir <- file.path("coverage_plots_difcover", species_name)
  setwd(working_dir)

  # Load all .DNAcopyout files
  files <- list.files(pattern = "\\.DNAcopyout$")
  together <- map_dfr(files, ~ read_tsv(.x, col_names = FALSE) %>%
    rename(scaf = X1, base = X2, stop = X3, windows = X4, log2_ratio = X5) %>%
    mutate(`bases spanned` = stop - base,
           scaf = gsub("^Chr", "", scaf)))

  # Load scaffold lengths
  scaffold_lengths <- read_tsv(ref_file, col_names = FALSE) %>%
    select(scaf = X1, LENGTH = X2) %>%
    mutate(scaf = gsub("^Chr", "", scaf)) %>%
    arrange(as.numeric(scaf)) %>%
    mutate(cum_start = lag(cumsum(LENGTH), default = 0))

  # Merge coverage data with scaffold lengths
  proportion <- together %>%
    left_join(scaffold_lengths, by = "scaf") %>%
    mutate(Bpos = cum_start + base,
           proportion = `bases spanned` / LENGTH)

  # Colors for this species
  colors <- color_map[[species_name]]
  if (is.null(colors)) stop("No colors defined for ", species_name)

  # Axis labels for Manhattan plot
  axis_set <- scaffold_lengths %>%
    mutate(center = cum_start + LENGTH / 2)

  # Scatter plot: log2 ratio vs proportion
  scatter_plot <- ggscatter(
    proportion,
    x = "log2_ratio", y = "proportion", color = "scaf", size = 2
  ) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    ggtitle(paste("DifCover scatter:", species_name, haplotype)) +
    theme_minimal()

  # Dot chart: chromosome-level % with significant coverage difference
  filtered_proportion <- proportion %>% filter(abs(log2_ratio) >= 0.7369656)
  chr_percentages <- filtered_proportion %>%
    group_by(scaf) %>%
    summarise(percent_diff_coverage = sum(`bases spanned`) / first(LENGTH),
              .groups = "drop")

  dot_chart <- ggdotchart(
    chr_percentages,
    x = "scaf",
    y = "percent_diff_coverage",
    color = colors[2],
    add = "segments",
    add.params = list(color = colors[1], size = 1),
    dot.size = 4
  ) +
    ggtitle(paste("Chromosomes with significantly different coverage:", species_name, haplotype)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Manhattan plot
  manhattan_plot <- ggplot(proportion, aes(x = Bpos, y = log2_ratio,
                                           color = as.factor(as.numeric(scaf) %% 2))) +
    geom_point(size = 1.2) +
    scale_color_manual(values = colors) +
    scale_x_continuous(labels = axis_set$scaf, breaks = axis_set$center) +
    ggtitle(paste("Manhattan plot:", species_name, haplotype)) +
    labs(x = "Chromosome", y = "Log2 ratio") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))

  # Save plots and data
  out_base <- file.path("..", "Results_plots", paste0(species_name, "_", haplotype, "_DifCover"))
  ggsave(paste0(out_base, ".pdf"), plot = ggarrange(scatter_plot, dot_chart, manhattan_plot, ncol = 1),
         width = 12, height = 16)
  ggsave(paste0(out_base, ".png"), plot = ggarrange(scatter_plot, dot_chart, manhattan_plot, ncol = 1),
         width = 12, height = 16, dpi = 300)

  write_tsv(filtered_proportion, paste0(out_base, "_significant_windows.tsv"))
  write_tsv(chr_percentages, paste0(out_base, "_chr_percent_diff_coverage.tsv"))
}

# Species and haplotypes to process
species_names <- c("RhinePhoxinus", "DanubeCsikii")
haplotypes <- c("Haplome_1", "Haplome_2")

walk2(
  rep(species_names, each = length(haplotypes)),
  rep(haplotypes, times = length(species_names)),
  plot_difcover_results
)
