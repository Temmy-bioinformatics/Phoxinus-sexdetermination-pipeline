library(tidyverse)
library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)

frq_file <- args[1]
indv_file <- args[2]
ldepth_mean_file <- args[3]
ldepth_file <- args[4]
imiss_file <- args[5]
lmiss_file <- args[6]
het_file <- args[7]
output_file <- args[8]

# 1. Variant Quality Check
variant_quality <- read_delim(ldepth_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)
p1 <- ggplot(variant_quality, aes(x = QUAL)) +
  geom_density(fill = "dodgerblue1", color = "black", alpha = 0.5) +
  labs(title = "Variant Quality Scores", x = "QUAL") +
  theme_minimal(base_size = 14) +
  coord_cartesian(xlim = c(0, 1000))

# 2. Mean Variant Depth
mean_depth <- read_delim(ldepth_mean_file, delim = "\t",
                         col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1, trim_ws = TRUE)
p2 <- ggplot(mean_depth, aes(x = mean_depth)) +
  geom_density(fill = "dodgerblue1", color = "black", alpha = 0.3) +
  labs(title = "Mean Depth Distribution", x = "Mean Depth") +
  theme_minimal(base_size = 14) +
  coord_cartesian(xlim = c(0, 50))

# 3. Variant Missingness
variant_missingness <- read_delim(lmiss_file, delim = "\t",
                                  col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1, trim_ws = TRUE)
p3 <- ggplot(variant_missingness, aes(x = fmiss)) +
  geom_density(fill = "dodgerblue1", color = "black", alpha = 0.3) +
  labs(title = "Variant Missingness", x = "Fraction Missing") +
  theme_minimal(base_size = 14)

# 4. Minor Allele Frequency (MAF)
allele_frequency <- read_delim(frq_file, delim = "\t",
                               col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1, trim_ws = TRUE)

# Calculate MAF assuming a1, a2 are allele counts; if frequencies, skip dividing by nchr
allele_frequency <- allele_frequency %>%
  mutate(maf = pmin(a1, a2) / nchr) %>%
  filter(!is.na(maf))

p4 <- ggplot(allele_frequency, aes(x = maf)) +
  geom_density(fill = "dodgerblue1", color = "black", alpha = 0.3) +
  labs(title = "Minor Allele Frequency (MAF)", x = "MAF") +
  theme_minimal(base_size = 14)

# 5. Mean Depth Per Individual
individual_depth <- read_delim(indv_file, delim = "\t",
                               col_names = c("ind", "nsites", "depth"), skip = 1, trim_ws = TRUE)
p5 <- ggplot(individual_depth, aes(x = depth)) +
  geom_histogram(fill = "dodgerblue1", color = "black", alpha = 0.3, bins = 30) +
  labs(title = "Mean Depth per Individual", x = "Depth") +
  theme_minimal(base_size = 14)

# 6. Missing Data Per Individual
individual_missingness <- read_delim(imiss_file, delim = "\t",
                                     col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1, trim_ws = TRUE)
p6 <- ggplot(individual_missingness, aes(x = fmiss)) +
  geom_histogram(fill = "dodgerblue1", color = "black", alpha = 0.3, bins = 30) +
  labs(title = "Missing Data per Individual", x = "Fraction Missing") +
  theme_minimal(base_size = 14)

# 7. Individual Heterozygosity
individual_heterozygosity <- read_delim(het_file, delim = "\t",
                                        col_names = c("ind", "ho", "he", "nsites", "f"), skip = 1, trim_ws = TRUE)
p7 <- ggplot(individual_heterozygosity, aes(x = f)) +
  geom_histogram(fill = "dodgerblue1", color = "black", alpha = 0.3, bins = 30) +
  labs(title = "Individual Heterozygosity (f)", x = "Heterozygosity (f)") +
  theme_minimal(base_size = 14)

# Combine plots
combined_plot <- (p1 | p2 | p3) / (p4 | p5 | p6) / p7

# Save to file
ggsave(output_file, combined_plot, width = 13.33, height = 7.5, dpi = 300)

print(combined_plot)
