# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3) {
  stop("Usage: Rscript plot_nuc_diversity.R <base_nuc_diversity_dir> <ref_file> <species>")
}

# Input directories and parameters
base_nuc_diversity_dir <- args[1]
ref_file <- args[2]
species <- args[3]

# Define color pairs for each species
color_map <- list(
  RhinePhoxinus = c("#9897A9", "#66A61E"),
  DanubeCsikii = c("#D95E02", "#9897A9")
)

# Use species-specific colors or default
cols <- if (species %in% names(color_map)) {
  color_map[[species]]
} else {
  c("#f8cfbeff", "#f04d29ff")  # Default colors
}

# Species directory path
species_dir <- file.path(base_nuc_diversity_dir, species)

# Load male and female nucleotide diversity data
mpi10_path <- file.path(species_dir, paste0(species, "_Male_10k_nucleotide_diversity.windowed.pi"))
fpi10_path <- file.path(species_dir, paste0(species, "_Female_10k_nucleotide_diversity.windowed.pi"))

if (!file.exists(mpi10_path) || !file.exists(fpi10_path)) {
  stop("Input nucleotide diversity files for male or female not found.")
}

mpi10 <- read_delim(mpi10_path, show_col_types = FALSE)
fpi10 <- read_delim(fpi10_path, show_col_types = FALSE)

# Join data and calculate differences
sex_pi10 <- full_join(mpi10, fpi10, by = c("CHROM", "BIN_START", "BIN_END")) %>%
  rename(
    scaf = CHROM,
    n_var_M = N_VARIANTS.x,
    pi_M = PI.x,
    n_var_F = N_VARIANTS.y,
    pi_F = PI.y
  ) %>%
  mutate(
    scaf = as.numeric(gsub("Chr", "", scaf)),
    pi_dif_M_minus_F = pi_M - pi_F,
    var_dif_M_minus_F = n_var_M - n_var_F,
    color_group = ifelse(scaf %% 2 == 0, "Even", "Odd")
  )

# Write differences to output file
write_tsv(sex_pi10, file.path(species_dir, paste0(species, "_NucDiv_10kb_sex_differences.txt")))

# Load reference genome info and calculate cumulative start positions
ref <- read.table(ref_file, header = FALSE, col.names = c("scaf", "LENGTH", "X1", "X2", "X3")) %>%
  mutate(scaf = as.numeric(gsub("Chr", "", scaf))) %>%
  arrange(scaf) %>%
  mutate(cum_start = cumsum(c(0, LENGTH[-n()])))

# Calculate chromosome midpoint for x-axis labels
axis.set <- ref %>%
  group_by(scaf) %>%
  summarize(center = (cum_start + (cum_start + LENGTH)) / 2, .groups = "drop")

# Merge cumulative start with diversity data
sex_pi10 <- left_join(sex_pi10, ref, by = "scaf") %>%
  arrange(scaf)

# Calculate genome-wide position for plotting
sex_pi10 <- sex_pi10 %>%
  mutate(Wpos = ((BIN_END - BIN_START) / 2) + BIN_START + cum_start)

# Plot variant count difference (male - female)
p_var <- ggplot(sex_pi10, aes(x = Wpos, y = var_dif_M_minus_F)) +
  geom_point(aes(color = color_group), size = 1) +
  scale_x_continuous(labels = axis.set$scaf, breaks = axis.set$center) +
  scale_color_manual(values = cols, guide = "none") +
  labs(x = "Chromosome", y = "Variants per 10kb (male - female)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(file.path(species_dir, paste0(species, "_var_per_10kb_window.tiff")), plot = p_var, width = 12, height = 12)

# Plot nucleotide diversity difference (pi) (male - female)
p_pi <- ggplot(sex_pi10, aes(x = Wpos, y = pi_dif_M_minus_F)) +
  geom_point(aes(color = color_group), size = 1) +
  scale_x_continuous(labels = axis.set$scaf, breaks = axis.set$center) +
  scale_color_manual(values = cols, guide = "none") +
  labs(x = "Chromosome", y = expression(pi ~ "(male mean - female mean)")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(file.path(species_dir, paste0(species, "_pi_10kb_window.tiff")), plot = p_pi, width = 12, height = 12)

cat("Nucleotide diversity analysis for species", species, "complete.\n")
