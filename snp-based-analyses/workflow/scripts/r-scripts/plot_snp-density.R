library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggpubr)
library(ggplot2)

# Define color pairs for each species
color_map <- list(
    SeinePhoxinus = c("#9897A9", "#ab9519")
)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3) {
  stop("Usage: Rscript plot_snp_density.R <base_dir> <ref_file> <species>")
}

base_dir <- args[1]
ref_file <- args[2]
species <- args[3]

# Get colors for species or default
colors <- color_map[[species]]
if (is.null(colors)) {
  colors <- c("skyblue", "orange") # default fallback colors
}

# Define output paths
output_dir <- file.path(base_dir, species, "plot")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

output_pdf <- file.path(output_dir, paste0(species, "_proportion_p001.pdf"))
output_txt <- file.path(output_dir, paste0(species, "_SNPdensity_SexFindR.txt"))
output_significant_txt <- file.path(output_dir, paste0(species, "_significant_SNPdensity_SexFindR.txt"))

# Function to load and prepare SNP density data for given sex
load_snp_density_data <- function(sex) {
  files <- list.files(path = file.path(base_dir, species), pattern = paste0(sex, ".*\\.snpden$"), full.names = TRUE)
  if (length(files) == 0) stop(paste("No SNP density files found for", sex, "in species", species))
  
  # Read the first file to create backbone
  backbone <- read_delim(files[1], delim = "\t", col_names = TRUE, show_col_types = FALSE) %>%
    mutate(CHROM = as.numeric(gsub("Chr", "", CHROM))) %>%
    unite(LOCATION, c("CHROM", "BIN_START"), sep = ":") %>%
    rename(!!sex := `VARIANTS/KB`) %>%
    select(-SNP_COUNT)
  
  # Merge with other files if multiple exist
  if (length(files) > 1) {
    for (file in files[-1]) {
      df <- read_delim(file, delim = "\t", col_names = TRUE, show_col_types = FALSE) %>%
        mutate(CHROM = as.numeric(gsub("Chr", "", CHROM))) %>%
        unite(LOCATION, c("CHROM", "BIN_START"), sep = ":") %>%
        rename(!!sex := `VARIANTS/KB`) %>%
        select(-SNP_COUNT)
      
      backbone <- full_join(backbone, df, by = "LOCATION")
    }
  }
  
  backbone %>%
    separate(LOCATION, into = c("scaf", "base"), sep = ":") %>%
    mutate(scaf = as.numeric(scaf), base = as.numeric(base))
}

# Load SNP density data for males and females
SNPdensity_males <- load_snp_density_data("Male")
SNPdensity_females <- load_snp_density_data("Female")

# Merge male and female data; replace NA with 0
SNPdensity <- full_join(SNPdensity_males, SNPdensity_females, by = c("scaf", "base")) %>%
  mutate(across(starts_with("Male"), ~replace_na(.x, 0)),
         across(starts_with("Female"), ~replace_na(.x, 0)))

# Calculate mean variant density for males and females and their difference
male_cols <- grep("Male", colnames(SNPdensity), value = TRUE)
female_cols <- grep("Female", colnames(SNPdensity), value = TRUE)

SNPdensity <- SNPdensity %>%
  mutate(
    mean_Males = rowMeans(across(all_of(male_cols)), na.rm = TRUE),
    mean_Females = rowMeans(across(all_of(female_cols)), na.rm = TRUE),
    mean_MvF_dif = mean_Males - mean_Females
  )

# Permutation test to compute p-values per position
set.seed(123)
perm_results <- replicate(1000, {
  permuted_cols <- sample(c(male_cols, female_cols))
  perm_males <- SNPdensity %>% select(all_of(permuted_cols[seq_along(male_cols)])) %>% rowMeans(na.rm = TRUE)
  perm_females <- SNPdensity %>% select(all_of(permuted_cols[(length(male_cols)+1):length(permuted_cols)])) %>% rowMeans(na.rm = TRUE)
  perm_males - perm_females
}, simplify = FALSE)

# Calculate empirical p-values
perm_with_true <- SNPdensity %>%
  mutate(
    Pvalue = sapply(seq_len(nrow(.)), function(i) {
      true_diff <- mean_MvF_dif[i]
      more_extreme <- sum(sapply(perm_results, function(perm) abs(perm[i]) >= abs(true_diff)))
      (more_extreme + 1) / (length(perm_results) + 1)
    })
  )

# Load reference genome and compute cumulative chromosome start positions
ref <- read.table(ref_file, header = FALSE, col.names = c("scaf", "LENGTH", "X1", "X2", "X3")) %>%
  mutate(scaf = as.numeric(gsub("Chr", "", scaf))) %>%
  arrange(scaf) %>%
  mutate(cum_start = cumsum(c(0, LENGTH[-n()])))

# Calculate proportion of significant SNPs (p <= 0.001) per scaffold
significant_snps <- perm_with_true %>%
  filter(Pvalue <= 0.001) %>%
  group_by(scaf) %>%
  summarize(n = n(), .groups = "drop") %>%
  left_join(ref, by = "scaf") %>%
  mutate(proportion = (n * 10000) / LENGTH)

# Save results
write_tsv(perm_with_true %>% select(scaf, base, mean_MvF_dif, Pvalue), output_txt)
write_tsv(significant_snps, output_significant_txt)

# Plot proportion of significant SNPs per scaffold
plot <- ggplot(significant_snps, aes(x = as.numeric(scaf), y = proportion)) +
  geom_col(aes(fill = as.factor(scaf %% 2)), show.legend = FALSE) +
  scale_fill_manual(values = colors) +
  labs(x = "Chromosome", y = "Proportion of significant SNPs (p ≤ 0.001)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(output_pdf, plot, width = 10, height = 6)

cat("SNP density analysis and plotting complete for species", species, "\n")
