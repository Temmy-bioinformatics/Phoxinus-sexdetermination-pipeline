# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

# Snakemake provides arguments for directories, files, and species
args <- commandArgs(trailingOnly = TRUE)

# Input directories and files
base_nuc_diversity_dir <- args[1]  # Base nucleotide diversity directory
ref_file <- args[2]  # Reference genome file
species <- args[3]  # Species name

# Define color pairs for each species
color_map <- list(
  "RhinePhoxinus" = c("#9897A9", "#66A61E"),
  "DanubeCsikii" = c("#D95E02", "#9897A9")
)

# Default color if species is not in the list
cols <- if (species %in% names(color_map)) {
  color_map[[species]]
} else {
  c("#f8cfbeff", "#f04d29ff")  # Default colors
}

# Set the working directory to the species' directory
species_dir <- file.path(base_nuc_diversity_dir, species)
setwd(species_dir)

# Load the data for males and females
mpi10 <- read_delim(paste0(species, "_Male_10k_nucleotide_diversity.windowed.pi"))
fpi10 <- read_delim(paste0(species, "_Female_10k_nucleotide_diversity.windowed.pi"))

# Merge male and female data, calculate differences, and prepare the data for plotting
sex_pi10 <- full_join(mpi10, fpi10, by = c("CHROM", "BIN_START", "BIN_END")) %>%
  rename(scaf = CHROM,  # Rename CHROM to scaf
         n_var_M = N_VARIANTS.x,
         pi_M = PI.x,
         n_var_F = N_VARIANTS.y,
         pi_F = PI.y) %>%
  mutate(
    scaf = as.numeric(gsub('Chr', '', scaf)),  # Remove "Chr" and convert to numeric
    pi_dif_M_minus_F = pi_M - pi_F,
    var_dif_M_minus_F = n_var_M - n_var_F,
    color_group = ifelse(scaf %% 2 == 0, "Even", "Odd")  # Add color_group column
  )

# Write the output file for sex differences in nucleotide diversity
write_tsv(sex_pi10, paste0(species, "_NucDiv_10kb_sex_differences.txt"))

# Load reference genome and calculate chromosome lengths
ref <- read.table(ref_file, header = FALSE, col.names = c("scaf", "LENGTH", "X1", "X2", "X3")) %>%
  mutate(scaf = as.numeric(gsub('Chr', '', scaf))) %>%
  arrange(scaf) %>%
  mutate(cum_start = cumsum(c(0, .data$LENGTH[-n()])))

# Midpoints for x axis
axis.set <- ref %>% 
  group_by(scaf) %>% 
  summarize(center = ((((cum_start + LENGTH) - cum_start) / 2) + cum_start))

# Merge the sex_pi10 data with reference data (for positions)
sex_pi10 <- left_join(sex_pi10, ref) %>%
  arrange(scaf)

# Calculate cumulative bases for each position in the genome
sex_pi10$Wpos <- ((sex_pi10$BIN_END - sex_pi10$BIN_START) / 2) + sex_pi10$BIN_START + sex_pi10$cum_start

##### SCATTER PLOT: Variant Difference (Male - Female)
# Plot variant counts difference
p_var <- sex_pi10 %>%
  ggplot(aes(x = Wpos, y = var_dif_M_minus_F)) +
  geom_point(aes(color = color_group), size = 1) +  # Use color_group for color
  labs(x = "Chromosome", y = "Variants per 10kb (male - female)") +
  scale_x_continuous(labels = axis.set$scaf, breaks = axis.set$center) +
  scale_color_manual(values = cols, guide = "none") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))

# Save the variant difference plot
ggsave(file.path(base_nuc_diversity_dir, species, paste0(species, "_var_per_10kb_window.tiff")), plot = p_var, width = 12, height = 12)

##### SCATTER PLOT: Pi Difference (Male - Female)
# Plot pi difference (Male - Female)
p_pi <- sex_pi10 %>%
  ggplot(aes(x = Wpos, y = pi_dif_M_minus_F)) +
  geom_point(aes(color = color_group), size = 1) +  # Use color_group for color
  labs(x = "Chromosome", y = expression(kb ~ pi ~ "(male mean - female mean)")) +
  scale_x_continuous(labels = axis.set$scaf, breaks = axis.set$center) +
  scale_color_manual(values = cols, guide = "none") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))

# Save the pi difference plot
ggsave(file.path(base_nuc_diversity_dir, species, paste0(species, "_pi_10kb_window.tiff")), plot = p_pi, width = 12, height = 12)

# Return success
cat("Nucleotide diversity analysis for species", species, "complete.\n")
