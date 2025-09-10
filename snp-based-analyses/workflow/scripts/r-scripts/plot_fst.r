# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

# Function to process each species directory for 10k windows
process_species_dir_10k <- function(base_dir, species_name, ref_file, color1, color2) {
  # Construct the directory path for the species
  species_dir <- file.path(base_dir, species_name)
  
  # Construct the file path for 10k windows Fst data
  fst_file <- file.path(species_dir, paste0(species_name, "_Fst_10k.windowed.weir.fst"))

  # Check if the file exists
  if (file.exists(fst_file)) {
    # Read and preprocess Fst data
    Fst <- read_tsv(fst_file, show_col_types = FALSE) %>%
      tidyr::replace_na(list(WEIGHTED_FST = 0)) %>%
      rename(scaf = CHROM, base = BIN_START) %>%
      filter(WEIGHTED_FST >= 0) %>%
      mutate(base = as.numeric(base),
             scaf = as.numeric(gsub("Chr", "", scaf)))

    # Load the .fai file to get chromosome lengths and preprocess
    ref <- read.table(ref_file, header = FALSE, col.names = c("scaf", "LENGTH", "X1", "X2", "X3")) %>%
      mutate(scaf = as.numeric(gsub("Chr", "", scaf)),
             LENGTH = as.numeric(LENGTH)) %>%
      arrange(scaf) %>%
      mutate(cum_start = cumsum(c(0, LENGTH[-n()]))) %>%
      mutate(color_group = factor(scaf %% 2))  # Alternating colors for scaffolds

    # Merge Fst with reference data
    Fst <- left_join(ref, Fst, by = "scaf") %>%
      arrange(scaf, base) %>%
      mutate(Bpos = cum_start + base)

    # Calculate cutoffs for top 5% and 1% of Fst values
    cutoff1 <- quantile(Fst$WEIGHTED_FST, 0.95, na.rm = TRUE)
    cutoff2 <- quantile(Fst$WEIGHTED_FST, 0.99, na.rm = TRUE)

    # Save cutoff values to a file
    threshold_file <- file.path(species_dir, paste0(species_name, "_Fst_10k_thresholds.txt"))
    write.table(
      data.frame(
        Threshold = c("Top 5%", "Top 1%"),
        Value = c(cutoff1, cutoff2)
      ),
      file = threshold_file,
      quote = FALSE,
      row.names = FALSE,
      sep = "\t"
    )

    # Create axis set for plotting
    axis.set <- Fst %>%
      group_by(scaf) %>%
      summarize(center = (max(Bpos) + min(Bpos)) / 2, .groups = "drop")

    # Define colors for the plot
    cols <- c("0" = color1, "1" = color2)

    # Create and save plot
    Fst_plot <- Fst %>%
      ggplot(aes(x = Bpos, y = WEIGHTED_FST, color = color_group)) +
      geom_point(size = 0.5) +
      scale_color_manual(values = cols, guide = "none") +
      scale_x_continuous(labels = axis.set$scaf, breaks = axis.set$center) +
      geom_hline(yintercept = cutoff2, linetype = "dotted", color = "black", size = 0.75) +
      geom_hline(yintercept = cutoff1, linetype = "dotted", color = "black", size = 1.5) +
      labs(
        title = paste("Fst for", species_name, "10k windows with 5% and 1% cutoff"),
        x = "CHROMOSOME",
        y = "WEIGHTED FST"
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )

    ggsave(file.path(species_dir, paste0(species_name, "_Fst_10k.tiff")), width = 12, height = 6, plot = Fst_plot)

  } else {
    cat("10k Fst file not found for species:", species_name, "\n")
  }
}

# Main function
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    stop("Usage: Rscript plot_10k_fst.R <base_dir> <ref_file> <species_name>")
  }
  base_dir <- args[1]
  ref_file <- args[2]
  species_name <- args[3]

  # Define color map
  color_map <- list(
    RhinePhoxinus = c("#9897A9", "#66A61E"),
    DanubeCsikii = c("#D95E02", "#9897A9")
  )

  if (!(species_name %in% names(color_map))) {
    stop("Species name not recognized.")
  }
  color1 <- color_map[[species_name]][1]
  color2 <- color_map[[species_name]][2]

  cat("Processing species:", species_name, "with colors:", color1, color2, "\n")

  process_species_dir_10k(base_dir, species_name, ref_file, color1, color2)
}

# Run main
main()
