# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Capture arguments for species and colors
base_dir <- args[1]
ref_file <- args[2]
species_name <- args[3]

# Define color pairs for each species
color_map <- list(
  "WeserMorella" = c("#da149b", "#9897A9"),
  "RhineCsikii" = c("#9897A9", "#654CFF"),
  "RhinePhoxinus" = c("#9897A9", "#66A61E"),
  "DanubeCsikii" = c("#D95E02", "#9897A9")
)

# Get the colors for the species, default to stop if species not recognized
if (species_name %in% names(color_map)) {
  color1 <- color_map[[species_name]][1]
  color2 <- color_map[[species_name]][2]
} else {
  stop("Species name not recognized.")
}

cat("Processing species:", species_name, "with colors:", color1, color2, "\n")

# Function to process each species directory
process_species_dir <- function(base_dir, species_name, ref_file, color1, color2) {
  # Construct the directory path for the species
  species_dir <- file.path(base_dir, species_name)
  
  # Construct the file path for Fst data
  fst_file <- file.path(species_dir, paste0(species_name, "_Fst.txt.weir.fst"))

  # Check if the file exists
  if (file.exists(fst_file)) {
    # Read and preprocess Fst data
    Fst_start <- read_tsv(fst_file, show_col_types = FALSE) %>%
      tidyr::replace_na(list(WEIGHTED_FST = 0)) %>%
      rename(scaf = CHROM, base = POS) %>%
      mutate(scaf = as.numeric(gsub('Chr', '', scaf)),
             base = as.numeric(base)) %>%
      filter(WEIR_AND_COCKERHAM_FST >= 0) %>%
      tidyr::replace_na(list(WEIR_AND_COCKERHAM_FST = 0))

    # Load the .fai file to get chromosome lengths and preprocess
    ref <- read.table(ref_file, header = FALSE, col.names = c("scaf", "LENGTH", "X1", "X2", "X3")) %>%
      mutate(scaf = as.numeric(gsub('Chr', '', scaf)),
             LENGTH = as.numeric(LENGTH)) %>% 
      arrange(scaf) %>%
      mutate(cum_start = cumsum(c(0, LENGTH[-n()]))) %>% 
      mutate(color_group = factor(scaf %% 2))  # Alternating colors for scaffolds

    # Merge Fst with reference data
    Fst <- left_join(ref, Fst_start, by = "scaf") %>%
      filter(WEIR_AND_COCKERHAM_FST >= 0) %>%
      arrange(scaf, base) %>%
      mutate(Bpos = cum_start + base)

    # Calculate cutoffs for top 5% and 1% of Fst values
    cutoff1 <- round(nrow(Fst) * 0.05)
    Fst_sort_cut1 <- Fst %>%
      arrange(desc(WEIR_AND_COCKERHAM_FST)) %>%
      slice_head(n = cutoff1)
    min_cutoff1 <- min(Fst_sort_cut1$WEIR_AND_COCKERHAM_FST)

    cutoff2 <- round(nrow(Fst) * 0.01)
    Fst_sort_cut2 <- Fst %>%
      arrange(desc(WEIR_AND_COCKERHAM_FST)) %>%
      slice_head(n = cutoff2)
    min_cutoff2 <- min(Fst_sort_cut2$WEIR_AND_COCKERHAM_FST)

    # Summary statistics
    quantiles <- quantile(Fst$WEIR_AND_COCKERHAM_FST, c(0.95, 0.99), na.rm = TRUE)
    mean_fst <- mean(Fst$WEIR_AND_COCKERHAM_FST, na.rm = TRUE)

    # Create axis set for plotting
    axis.set <- Fst %>%
      group_by(scaf) %>%
      summarize(center = (max(Bpos) + min(Bpos)) / 2, .groups = 'drop')

    # Define colors for the plot
    cols <- c("0" = color1, "1" = color2)  # Use passed colors

    # Create and save plots
    Fst_plot <- Fst %>%
      ggplot(aes(x = Bpos, y = WEIR_AND_COCKERHAM_FST, color = color_group)) +
      geom_point(size = 0.5) +  # Points for each data point
      scale_color_manual(values = cols, guide = "none") +  # Use passed colors
      scale_x_continuous(labels = axis.set$scaf, breaks = axis.set$center) +  # Chromosome labels
      geom_hline(yintercept = min_cutoff2, linetype = "dotted", color = "black", size = 0.75) +  # 1% cutoff line
      geom_hline(yintercept = min_cutoff1, linetype = "dotted", color = "black", size = 1.5) +  # 5% cutoff line
      labs(
        title = paste("Fst for", species_name, "with 5% and 1% cutoff as dotted lines"),
        x = "CHROMOSOME",
        y = "WEIR AND COCKERHAM FST"
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate chromosome labels for clarity
        panel.grid.major.x = element_blank(),  # Remove vertical grid lines for cleaner look
        panel.grid.minor.x = element_blank()
      )

    # Save the plot
    ggsave(file.path(species_dir, paste0(species_name, "_Fst.tiff")), width = 12, height = 6, plot = Fst_plot)

    # Prepare input file for combined analysis
    Fst_sort_cut <- Fst %>%
      arrange(desc(WEIR_AND_COCKERHAM_FST)) %>%
      slice_head(n = cutoff1) %>%
      select(scaf, base)
    write_tsv(Fst_sort_cut, path = file.path(species_dir, paste0(species_name, "_all_Fst_sort_cut.txt")), col_names = FALSE)

  } else {
    cat("Fst file not found for species:", species_name, "\n")
  }
}

# Main function to process one species
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  base_dir <- args[1]
  ref_file <- args[2]
  species_name <- args[3]  # Now we process one species at a time

  # Call the function to process the species
  process_species_dir(base_dir, species_name, ref_file, color1, color2)
}

# Run the main function
main()
