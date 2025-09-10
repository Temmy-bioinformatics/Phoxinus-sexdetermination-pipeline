# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# Define color pairs for each species
color_map <- list(
  "WeserMorella" = c("#da149b", "#9897A9"),
  "RhineCsikii" = c("#9897A9", "#654CFF"),
  "RhinePhoxinus" = c("#9897A9", "#66A61E"),
  "DanubeCsikii" = c("#D95E02", "#9897A9")
)


# Function to process each species directory
process_species_dir <- function(base_dir, species_name, ref_file) {
  # Construct the species-specific directory paths
  species_dir <- file.path(base_dir, species_name)
  output_dir <- file.path(species_dir, "output")
  
  # Print species information
  cat("Processing species:", species_name, "in directory:", output_dir, "\n")

  # Define colors based on species
  if (species_name %in% names(color_map)) {
    color1 <- color_map[[species_name]][1]
    color2 <- color_map[[species_name]][2]
  } else {
  }

  # Debug: Show selected colors for the species
  cat("Using colors:", color1, color2, "\n")

  # Construct the file path for GWAS data inside the output subdirectory
  gwas_file <- file.path(output_dir, paste0(species_name, "_GEMMA.assoc.txt"))

  # Check if the GWAS file exists
  if (!file.exists(gwas_file)) {
    cat("GWAS file not found for species:", species_name, "\n")
    return()
  }

  # Read in GWAS data
  gwas0 <- read_tsv(gwas_file, show_col_types = FALSE)

  # Debug: Check initial NA values in key columns
  cat("Initial NAs in GWAS data (ps):", sum(is.na(gwas0$ps)), "\n")
  cat("Initial NAs in GWAS data (p_lrt):", sum(is.na(gwas0$p_lrt)), "\n")

  # Rename the 'chr' column to 'scaf'
  gwas0 <- gwas0 %>%
    rename(scaf = chr)

  # Replace -nan with NA in the 'ps' and 'p_lrt' columns
  gwas1 <- gwas0 %>%
    mutate(
      ps = ifelse(is.nan(ps), NA, ps),
      p_lrt = ifelse(is.nan(p_lrt), NA, p_lrt)
    )

  # Debug: Check NAs after replacement
  cat("NAs after replacement in GWAS data (ps):", sum(is.na(gwas1$ps)), "\n")
  cat("NAs after replacement in GWAS data (p_lrt):", sum(is.na(gwas1$p_lrt)), "\n")

  # Select necessary columns and ensure they are numeric
  gwas1 <- gwas1 %>%
    select(scaf, ps, n_mis, n_obs, allele0, allele1, af, p_lrt) %>%
    mutate(
      scaf = as.numeric(scaf),  # Convert 'scaf' to character for consistent type
      ps = as.numeric(ps),
      p_lrt = as.numeric(p_lrt)
    ) %>%
    filter(!is.na(ps) & !is.na(p_lrt))  # Ensure ps and p_lrt are numeric

  # Debug: Check NAs after filtering
  cat("NAs after filtering in GWAS data (ps):", sum(is.na(gwas1$ps)), "\n")
  cat("NAs after filtering in GWAS data (p_lrt):", sum(is.na(gwas1$p_lrt)), "\n")

  # Find top 5% and 1% of p-values using the cutoff approach
  cutoff1 <- round(nrow(gwas1) * 0.05)  # Top 5%
  gwas_sort_cut1 <- gwas1 %>% arrange(p_lrt) %>% slice_head(n = cutoff1)
  perc5 <- max(gwas_sort_cut1$p_lrt)
  cat("Top 5% p-value cutoff:", perc5, "\n")

  cutoff2 <- round(nrow(gwas1) * 0.01)  # Top 1%
  gwas_sort_cut2 <- gwas1 %>% arrange(p_lrt) %>% slice_head(n = cutoff2)
  perc1 <- max(gwas_sort_cut2$p_lrt)
  cat("Top 1% p-value cutoff:", perc1, "\n")

  # Convert the p-values to the -log10 scale
  gwas1 <- gwas1 %>%
    mutate(p_lrt = -log10(p_lrt))

  # Load the reference file (.fai) to get chromosome lengths and preprocess
  ref <- read.table(ref_file, header = FALSE, col.names = c("scaf", "LENGTH", "X1", "X2", "X3")) %>%
      mutate(scaf = as.numeric(gsub('Chr', '', scaf)),
             LENGTH = as.numeric(LENGTH)) %>% 
      arrange(scaf) %>%
      mutate(cum_start = cumsum(c(0, LENGTH[-n()]))) %>% 
      mutate(color_group = factor(scaf %% 2))  # Alternating colors for scaffolds


 #make numeric
  gwas1$scaf <- as.numeric(gwas1$scaf)
  
  # Perform the join of the GWAS data with the reference data (including color_group)

  gwas <- left_join(gwas1, ref, by = "scaf") %>%
    arrange(scaf) %>%
    mutate(Bpos = cum_start + ps)  # Calculate Bpos for each position

  # Create axis set for plotting (center positions of each scaffold)
  axis_set <- gwas %>%
    group_by(scaf) %>%
    reframe(center = (max(Bpos) + min(Bpos)) / 2)

  # Debug: Check NAs after join
  cat("NAs after join in GWAS data (Bpos):", sum(is.na(gwas$Bpos)), "\n")
  # Debug: Check NAs after calculating Bpos
  cat("NAs after calculating Bpos in GWAS data:", sum(is.na(gwas$Bpos)), "\n")

  # Define colors for the plot (using species-specific colors)
  cols <- c("0" = color1, "1" = color2)

  # Create the plot
  gwas_plot <- ggplot(gwas, aes(x = Bpos, y = p_lrt, color = color_group)) +
    geom_point(size = 0.5) +  # Plot points for each data
    scale_color_manual(values = cols, guide = "none") +  # Use defined colors for alternating groups
    scale_x_continuous(labels = axis_set$scaf, breaks = axis_set$center) +  # X-axis breaks based on chromosome centers
    labs(x = "CHROMOSOME", y = "-log10 p-value") +  # Axis labels
    ggtitle(paste("GWAS outliers for", species_name, "with 5% and 1% cutoff as dotted lines")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))  # Rotate x-axis labels for clarity

  # Save the plot as both TIFF and PDF
  ggsave(file.path(output_dir, paste0(species_name, "_GWAS.tiff")), width = 12, height = 6, plot = gwas_plot)

  # Write the complete GWAS data to a file
  write_tsv(gwas, file = file.path(output_dir, paste0(species_name, "_gwas_plot_values.txt")), col_names = TRUE)

  # Prepare input file for combined analysis (top 5% of p-values)
  GWAS_out_sort_cut <- gwas0 %>%
    arrange(p_lrt) %>%
    slice_head(n = cutoff1) %>%
    select(scaf, ps)

  write_tsv(GWAS_out_sort_cut, file = file.path(output_dir, paste0(species_name, "_GWAS_out_sort_cut.txt")), col_names = FALSE)
}

# Main function to process one species
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  base_dir <- args[1]
  ref_file <- args[2]
  species_name <- args[3]  # Now we process one species at a time

  # Process the species directory
  process_species_dir(base_dir, species_name, ref_file)
}

# Run the main function
main()
