# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# Clear the workspace
rm(list = ls())

# Set the base directory where all subdirectories (species directories) are located
base_dir <- "/home/toriowo/SNAKEMAKE/SexDetermination_snakemake/results/SNP_Gwas/GWAS"

# Load the .fai file to get chromosome lengths and rename the 'CHR' column to 'scaf'
ref_file <- "/home/toriowo/Sex_Determination/Results/Final_Phoxy_1_rearranged.fasta.fai"
ref <- read.table(ref_file, header = FALSE, col.names = c("scaf", "LENGTH", "X1", "X2", "X3"))

# Calculate cumulative chromosome start positions
ref <- ref %>%
  mutate(scaf = as.numeric(gsub('Chr', '', scaf))) %>%
  arrange(scaf) %>%
  mutate(cum_start = cumsum(c(0, LENGTH[-nrow(ref)])))

# Get a list of all subdirectories (each representing a species)
species_dirs <- list.dirs(base_dir, recursive = FALSE)

# Loop through each species directory
for (species_dir in species_dirs) {
  # Define the output subdirectory path
  output_dir <- file.path(species_dir, "output")
  
  # Extract the species name from the directory path
  species_name <- basename(species_dir)
  cat("Processing species:", species_name, "in directory:", output_dir, "\n")

  # Construct the file path for GWAS data inside the output subdirectory
  gwas_file <- file.path(output_dir, paste0(species_name, "_GEMMA.assoc.txt"))

  # Check if the GWAS file exists
  if (file.exists(gwas_file)) {
    # Read in GWAS data
    gwas1 <- read_tsv(gwas_file)
    
    # Rename 'chr' to 'scaf'
    gwas2 <- gwas1 %>%
      rename(scaf = chr)

    # Select necessary columns
    gwas3 <- gwas2 %>%
      select(scaf, ps, n_mis, n_obs, allele0, allele1, af, p_lrt)

    # Find top 5% and 1% of p-values
    cutoff1 <- round(nrow(gwas3) * 0.05)
    gwas_sort_cut1 <- gwas3 %>%
      arrange(p_lrt) %>%
      slice_head(n = cutoff1)

    cutoff2 <- round(nrow(gwas3) * 0.01)
    gwas_sort_cut2 <- gwas3 %>%
      arrange(p_lrt) %>%
      slice_head(n = cutoff2)

    # Convert the p-values to the -log10 scale
    gwas4 <- gwas3 %>%
      mutate(p_lrt = -log10(p_lrt))

    # Ensure both columns are of type character
    gwas4 <- gwas4 %>%
      mutate(scaf = as.character(scaf))

    ref <- ref %>%
      mutate(scaf = as.character(scaf))

    # Perform the join
    gwas <- left_join(gwas4, ref, by = "scaf") %>%
      arrange(scaf)

    # For x-axis, use cumulative bases for each position in the genome
    gwas <- gwas %>%
      mutate(Bpos = cum_start + ps)

    # Set x-axis breaks and labels
    axis_set <- gwas %>%
      group_by(scaf) %>%
      summarize(center = (max(Bpos) + min(Bpos)) / 2, .groups = 'drop')

    # Define colors for the plot
    cols <- c("n" = "#d7edfcee", "y" = "#219af1ee")

    # #Uncomment and modify this section to create and save plots as needed
    # gwas_plot <- gwas %>%
    #   ggplot(aes(x = Bpos, y = p_lrt, colour = fill)) +
    #   geom_point(size = 1) +
    #   scale_color_manual(values = cols, guide = "none") +
    #   scale_x_continuous(labels = axis_set$scaf, breaks = axis_set$center) +
    #   labs(x = "CHROMOSOME", y = "-log10 p-value") +
    #   geom_hline(yintercept = -log10(0.1264462), linetype = "dotted", color = "black", size = 1.5) +
    #   geom_hline(yintercept = -log10(0.03907123), linetype = "dotted", color = "black", size = 0.75) +
    #   ggtitle(paste("GWAS outliers for", species_name, "with 5% and 1% cutoff as dotted lines")) +
    #   theme_bw() +
    #   theme(axis.text.x = element_text(angle = 90))

    # ggsave(file.path(output_dir, paste0(species_name, "_GWAS.tiff")), plot = gwas_plot)
    # ggsave(file.path(output_dir, paste0(species_name, "_GWAS.pdf")), plot = gwas_plot)

    # Prepare input file for combined analysis
    cutoff <- round(nrow(gwas3) * 0.05)
    GWAS_out_sort_cut <- gwas3 %>%
      arrange(p_lrt) %>%
      slice_head(n = cutoff) %>%
      select(scaf, ps)

    write_tsv(GWAS_out_sort_cut, path = file.path(output_dir, paste0(species_name, "_GWAS_out_sort_cut.txt")), col_names = FALSE)

  } else {
    cat("GWAS file not found for species:", species_name, "in directory:", output_dir, "\n")
  }
}

cat("Processing complete.\n")
