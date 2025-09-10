library(dplyr) 
library(readr)
library(stringr)
library(tidyr)
library(ggpubr)
library(ggplot2)

# Define color pairs for each species
color_map <- list(
  "WeserMorella" = c("#da149b", "#9897A9"),
  "RhineCsikii" = c("#9897A9", "#654CFF"),
  "RhinePhoxinus" = c("#9897A9", "#66A61E"),
  "DanubeCsikii" = c("#D95E02", "#9897A9")
)

# Command-line arguments for file paths and species name
args <- commandArgs(trailingOnly = TRUE)
base_dir <- args[1]
ref_file <- args[2]
species <- args[3]

# Get color pair for the species or use default colors
colors <- color_map[[species]]
if (is.null(colors)) {
  colors <- c("skyblue", "orange") # Default colors
}

# Derived output file paths
output_dir <- file.path(base_dir, species, "plot")
output_pdf <- file.path(output_dir, paste0(species, "_proportion_p001.pdf"))
output_txt <- file.path(output_dir, paste0(species, "_SNPdensity_SexFindR.txt"))
output_significant_txt <- file.path(output_dir, paste0(species, "_significant_SNPdensity_SexFindR.txt"))

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Function to load SNP density data
load_snp_density_data <- function(sex) {
  myFiles <- list.files(path = file.path(base_dir, species), pattern = paste0(sex, ".*\\.snpden$"), full.names = TRUE)
  if (length(myFiles) == 0) stop(paste("No files found for", sex, "in", species))
  
  backbone <- read_delim(file = myFiles[1], delim = "\t", col_names = TRUE)
  column_ID <- str_extract(myFiles[1], "Male|Female")
  
  backbone.upgrade <- backbone %>%
    mutate(CHROM = as.numeric(gsub("Chr", "", CHROM))) %>% 
    unite(LOCATION, c("CHROM", "BIN_START"), sep = ":") %>%
    rename(!!column_ID := `VARIANTS/KB`) %>%
    select(-SNP_COUNT)
  
  for (file in myFiles[-1]) {
    file_data <- read_delim(file, delim = "\t", col_names = TRUE)
    column_ID <- str_extract(file, "Male|Female")
    
    file.upgrade <- file_data %>%
      mutate(CHROM = as.numeric(gsub("Chr", "", CHROM))) %>%
      unite(LOCATION, c("CHROM", "BIN_START"), sep = ":") %>%
      rename(!!column_ID := `VARIANTS/KB`) %>%
      select(-SNP_COUNT)
    
    backbone.upgrade <- full_join(backbone.upgrade, file.upgrade, by = "LOCATION")
  }
  
  backbone.upgrade %>%
    separate(LOCATION, into = c("scaf", "base"), sep = ":") %>%
    mutate(scaf = as.numeric(scaf),  # Ensure scaf is numeric
           base = as.numeric(base))
}

# Load data for males and females
SNPdensity.males <- load_snp_density_data("Male")
SNPdensity.females <- load_snp_density_data("Female")

# Merge data and handle missing values
SNPdensity <- full_join(SNPdensity.males, SNPdensity.females) %>%
  mutate(across(starts_with("Male"), ~ replace_na(.x, 0)),
         across(starts_with("Female"), ~ replace_na(.x, 0)))

# Calculate mean values for males and females
male_cols <- grep("Male", names(SNPdensity), value = TRUE)
female_cols <- grep("Female", names(SNPdensity), value = TRUE)

SNPdensity <- SNPdensity %>%
  mutate(mean_Males = rowMeans(across(all_of(male_cols)), na.rm = TRUE),
         mean_Females = rowMeans(across(all_of(female_cols)), na.rm = TRUE),
         mean_MvF_dif = mean_Males - mean_Females)

# Permutation test setup
set.seed(123)
perm_results <- replicate(1000, {
  perm_indices <- sample(c(male_cols, female_cols))
  perm_males <- SNPdensity %>% select(all_of(perm_indices[1:length(male_cols)])) %>%
    rowMeans(na.rm = TRUE)
  perm_females <- SNPdensity %>% select(all_of(perm_indices[(length(male_cols) + 1):length(perm_indices)])) %>%
    rowMeans(na.rm = TRUE)
  perm_males - perm_females
}, simplify = FALSE)

# Calculate p-values for the true data and convert scaf column to numeric
perm_with_true <- SNPdensity %>%
  mutate(
    Pvalue = sapply(1:nrow(.), function(i) {
      true_dif <- .$mean_MvF_dif[i]
      more_extreme <- sum(sapply(perm_results, function(perm) abs(perm[i]) >= abs(true_dif)))
      (more_extreme + 1) / (length(perm_results) + 1)
    }),
    scaf = as.numeric(scaf)  # Convert scaf column to numeric
  )

# Extract significant SNPs with P-value <= 0.001
significant_snps <- perm_with_true %>% filter(Pvalue <= 0.001)

# Save full results and significant SNPs to file
write_tsv(perm_with_true %>% select(scaf, base, mean_MvF_dif, Pvalue), output_txt)

# Load reference genome and calculate chromosome lengths
ref <- read.table(ref_file, header = FALSE, col.names = c("scaf", "LENGTH", "X1", "X2", "X3")) %>%
  mutate(scaf = as.numeric(gsub('Chr', '', scaf))) %>%
  arrange(scaf) %>%
  mutate(cum_start = cumsum(c(0, .data$LENGTH[-n()])))

# Add proportions of significant SNPs per scaffold
significant_snps <- significant_snps %>%
  group_by(scaf) %>%
  summarize(n = n(), .groups = "drop") %>%
  left_join(ref, by = "scaf") %>%
  mutate(proportion = (n * 10000) / LENGTH)

# Save significant SNPs to a separate file
write_tsv(significant_snps, output_significant_txt)

# Create and save the plot
plot <- ggplot(significant_snps, aes(x = as.numeric(scaf), y = proportion)) +
  geom_col(aes(fill = as.factor(as.numeric(scaf) %% 2)), show.legend = FALSE) +
  scale_fill_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Chromosome", y = "Proportion")

ggsave(output_pdf, plot)
