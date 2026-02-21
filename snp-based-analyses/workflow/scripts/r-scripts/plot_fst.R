# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

# Main function
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Check for correct number of arguments
  if (length(args) != 5) {
    stop("Usage: Rscript plot_fst.R <fst_file> <ref_file> <species_name> <output_plot> <output_thresholds>")
  }
  
  # Parse arguments
  fst_file <- args[1]
  ref_file <- args[2]
  species_name <- args[3]
  output_plot <- args[4]
  output_thresholds <- args[5]

  # Define color map for different species
  color_map <- list(
    SeinePhoxinus = c("#9897A9", "#ab9519")
    # Add more species here as needed:
    # SpeciesName = c("color1", "color2")
  )

  # Validate species name
  if (!(species_name %in% names(color_map))) {
    stop("Species name '", species_name, "' not recognized. Available species: ", 
         paste(names(color_map), collapse = ", "))
  }
  
  # Get colors for this species
  color1 <- color_map[[species_name]][1]
  color2 <- color_map[[species_name]][2]

  # Log processing information
  cat("================================================\n")
  cat("Processing FST plot for:", species_name, "\n")
  cat("================================================\n")
  cat("Input FST file:       ", fst_file, "\n")
  cat("Reference file:       ", ref_file, "\n")
  cat("Output plot:          ", output_plot, "\n")
  cat("Output thresholds:    ", output_thresholds, "\n")
  cat("Colors:               ", color1, ",", color2, "\n")
  cat("================================================\n\n")

  # Check if input files exist
  if (!file.exists(fst_file)) {
    stop("FST file not found: ", fst_file)
  }
  
  if (!file.exists(ref_file)) {
    stop("Reference file not found: ", ref_file)
  }

  # Read and preprocess Fst data
  cat("Reading FST data...\n")
  Fst <- read_tsv(fst_file, show_col_types = FALSE) %>%
    tidyr::replace_na(list(WEIGHTED_FST = 0)) %>%
    rename(scaf = CHROM, base = BIN_START) %>%
    filter(WEIGHTED_FST >= 0) %>%
    mutate(
      base = as.numeric(base),
      scaf = as.numeric(gsub("Chr", "", scaf))
    )
  
  cat("  - Loaded", nrow(Fst), "FST windows\n")

  # Load the .fai file to get chromosome lengths and preprocess
  cat("Reading reference genome index...\n")
  ref <- read.table(ref_file, header = FALSE, 
                    col.names = c("scaf", "LENGTH", "X1", "X2", "X3")) %>%
    mutate(
      scaf = as.numeric(gsub("Chr", "", scaf)),
      LENGTH = as.numeric(LENGTH)
    ) %>%
    arrange(scaf) %>%
    mutate(
      cum_start = cumsum(c(0, LENGTH[-n()])),
      color_group = factor(scaf %% 2)  # Alternating colors for scaffolds
    )
  
  cat("  - Found", nrow(ref), "chromosomes\n")

  # Merge Fst with reference data
  cat("Merging FST data with reference positions...\n")
  Fst <- left_join(ref, Fst, by = "scaf") %>%
    arrange(scaf, base) %>%
    mutate(Bpos = cum_start + base)

  # Calculate cutoffs for top 5% and 1% of Fst values
  cat("Calculating FST thresholds...\n")
  cutoff1 <- quantile(Fst$WEIGHTED_FST, 0.95, na.rm = TRUE)
  cutoff2 <- quantile(Fst$WEIGHTED_FST, 0.99, na.rm = TRUE)
  
  cat("  - Top 5% threshold: ", round(cutoff1, 4), "\n")
  cat("  - Top 1% threshold: ", round(cutoff2, 4), "\n")

  # Save cutoff values to a file
  cat("Writing threshold file...\n")
  write.table(
    data.frame(
      Threshold = c("Top 5%", "Top 1%"),
      Value = c(cutoff1, cutoff2)
    ),
    file = output_thresholds,
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )
  cat("  ✓ Thresholds saved to:", output_thresholds, "\n")

  # Create axis set for plotting
  axis.set <- Fst %>%
    group_by(scaf) %>%
    summarize(center = (max(Bpos) + min(Bpos)) / 2, .groups = "drop")

  # Define colors for the plot
  cols <- c("0" = color1, "1" = color2)

  # Create plot
  cat("Generating plot...\n")
  Fst_plot <- Fst %>%
    ggplot(aes(x = Bpos, y = WEIGHTED_FST, color = color_group)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = cols, guide = "none") +
    scale_x_continuous(labels = axis.set$scaf, breaks = axis.set$center) +
    geom_hline(yintercept = cutoff2, linetype = "dotted", color = "black", linewidth = 0.75) +
    geom_hline(yintercept = cutoff1, linetype = "dotted", color = "black", linewidth = 1.5) +
    labs(
      title = paste("Fst for", species_name, "- 10k windows with 5% and 1% cutoffs"),
      x = "CHROMOSOME",
      y = "WEIGHTED FST"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )

  # Save plot
  ggsave(output_plot, width = 12, height = 6, plot = Fst_plot, dpi = 300)
  cat("  ✓ Plot saved to:", output_plot, "\n")
  
  cat("\n================================================\n")
  cat("FST plotting completed successfully!\n")
  cat("================================================\n")
}

# Run main
main()