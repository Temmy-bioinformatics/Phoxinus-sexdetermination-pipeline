# Load required libraries
library(tidyverse)
library(patchwork)
library(ggthemes)

# Clear workspace
rm(list = ls())

# ------------------------
# USER SETTINGS
# ------------------------

# Set species name (used to load relevant files)
species <- "DanubeCsikii"

# Set project directories (modify as needed)
project_dir <- "/PROJECTS/Phoxinus_SexDetermination_Publication"
results_dir <- file.path(project_dir, "Randomised_subsampling_results", species)
ref_file <- "resources/refgenome/Final_Phoxy_1_rearranged.fasta.fai"

# Set output directory for plots and results
output_dir <- file.path(project_dir, "Results_plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------
# HELPER FUNCTIONS
# ------------------------

# Function to add cumulative start position by chromosome (not used below but useful)
add_cumulative_start <- function(df) {
  df %>%
    arrange(scaf) %>%
    mutate(start = cumsum(Length) - Length)
}

# ------------------------
# LOAD AND PROCESS DATA
# ------------------------

# SNP density data
snp_file <- file.path(results_dir, paste0(species, "_SNPdensity_SexFindR.txt"))
SNP <- read_delim(snp_file) %>%
  select(scaf, base, mean_MvF_dif, Pvalue) %>%
  mutate(
    scaf = str_remove(scaf, "Chr"), # remove 'Chr' prefix
    base = base + 10000
  ) %>%
  filter(Pvalue <= 0.05)

# Rank SNP density windows
SNP_filter_arrange <- SNP %>%
  arrange(Pvalue, desc(abs(mean_MvF_dif))) %>%
  mutate(SNPdensity_rank = row_number(),
         scaf = as.numeric(scaf))

# GWAS window count
gwas_file <- file.path(results_dir, paste0(species, "_GWAS_window_count.txt"))
GWAS_window_count <- read_delim(gwas_file) %>%
  arrange(desc(count)) %>%
  mutate(GWAScount_rank = row_number(),
         scaf = as.numeric(str_remove(scaf, "Chr")))

# Fst window count
fst_file <- file.path(results_dir, paste0(species, "_Fst_window_count.txt"))
Fst_window_count <- read_delim(fst_file) %>%
  arrange(desc(count)) %>%
  mutate(Fstcount_rank = row_number())

# Nucleotide diversity (Pi) data
pi_file <- file.path(results_dir, paste0(species, "_NucDiv_10kb_sex_differences.txt"))
pi <- read_delim(pi_file) %>%
  select(scaf, BIN_START, pi_dif_M_minus_F) %>%
  arrange(desc(abs(pi_dif_M_minus_F))) %>%
  mutate(pi_rank = row_number(),
         base = BIN_START + 9999,
         scaf = as.numeric(scaf))

# Load reference chromosome lengths and calculate cumulative starts
ref <- read.table(ref_file, header = FALSE, col.names = c("scaf", "LENGTH", "X1", "X2", "X3")) %>%
  mutate(scaf = as.numeric(str_remove(scaf, "Chr"))) %>%
  arrange(scaf) %>%
  mutate(cum_start = cumsum(c(0, LENGTH[-length(LENGTH)])))

# Define colors for plotting chromosomes (choose one set)
color_pairs <- c("#9897A9", "#66A61E") # example: RhinePhoxinus colors

# Assign color groups alternating by chromosome
ref <- ref %>%
  mutate(color_group = color_pairs[(scaf %% 2) + 1])

# Midpoints for x-axis labeling
axis_set <- ref %>%
  group_by(scaf) %>%
  summarize(center = cum_start + LENGTH / 2)

# ------------------------
# JOIN DATASETS
# ------------------------

Fst_window_count_chr <- left_join(Fst_window_count, ref, by = "scaf")
GWAS_window_count_chr <- left_join(GWAS_window_count, ref, by = "scaf")
SNP_filter_arrange_chr <- left_join(SNP_filter_arrange, ref, by = "scaf")

SNP_Fst <- full_join(SNP_filter_arrange_chr, Fst_window_count_chr, by = c("scaf", "base"))
SNP_Fst_GWAS <- left_join(SNP_Fst, GWAS_window_count, by = c("scaf", "base"))

all <- full_join(SNP_Fst_GWAS, pi, by = c("scaf", "base")) %>%
  select(-count.x, -count.y) %>%
  mutate(total_rank = GWAScount_rank + SNPdensity_rank + Fstcount_rank + pi_rank)

# ------------------------
# SPEARMAN CORRELATION OF RANKINGS
# ------------------------

correlations <- c(
  "SNP_density vs GWAS_count" = cor(all$SNPdensity_rank, all$GWAScount_rank, method = "spearman", use = "complete.obs"),
  "SNP_density vs Fst_count" = cor(all$SNPdensity_rank, all$Fstcount_rank, method = "spearman", use = "complete.obs"),
  "SNP_density vs pi" = cor(all$SNPdensity_rank, all$pi_rank, method = "spearman", use = "complete.obs"),
  "GWAS_count vs Fst_count" = cor(all$GWAScount_rank, all$Fstcount_rank, method = "spearman", use = "complete.obs"),
  "GWAS_count vs pi" = cor(all$GWAScount_rank, all$pi_rank, method = "spearman", use = "complete.obs"),
  "Fst_count vs pi" = cor(all$Fstcount_rank, all$pi_rank, method = "spearman", use = "complete.obs")
)

# Write sorted correlations to file
sink(file.path(project_dir, paste0(species, "_correlation_rankings.txt")))
cat("Ranked Spearman's Correlations (from strongest to weakest):\n\n")
sorted_correlations <- sort(abs(correlations), decreasing = TRUE)
for (name in names(sorted_correlations)) {
  cat(name, " - rho = ", correlations[name], "\n")
}
sink()

# ------------------------
# FILTER TOP WINDOWS BY RANK (top 100)
# ------------------------

top_all <- all %>%
  filter(GWAScount_rank <= 100,
         Fstcount_rank <= 100,
         pi_rank <= 100,
         SNPdensity_rank <= 100)

# Example 3-way filters (can add more as needed)
top_rank_SFP <- all %>% filter(SNPdensity_rank <= 100, Fstcount_rank <= 100, pi_rank <= 100)
top_rank_GFP <- all %>% filter(GWAScount_rank <= 100, Fstcount_rank <= 100, pi_rank <= 100)

# Save filtered datasets
write_tsv(top_all %>% select(-total_rank), file.path(output_dir, "candidate_windows_SNP_GWAS_Fst_pi.txt"))
write_tsv(top_rank_GFP %>% select(-total_rank), file.path(output_dir, "candidate_windows_GWAS_Pi_Fst.txt"))
write_tsv(top_rank_SFP %>% select(-total_rank), file.path(output_dir, "candidate_windows_SNP_Fst_Pi.txt"))

# ------------------------
# PREPARE WINDOW POSITIONS FOR PLOTTING
# ------------------------

add_Wpos <- function(df) {
  df %>%
    mutate(Wpos = cum_start + base)
}

Fst_window_count_chr <- add_Wpos(Fst_window_count_chr)
GWAS_window_count_chr <- add_Wpos(GWAS_window_count_chr)
SNP_filter_arrange_chr <- add_Wpos(SNP_filter_arrange_chr)
pi_chr <- add_Wpos(left_join(pi, ref, by = "scaf"))

top_points <- add_Wpos(top_all)
top_window <- top_points$Wpos

# ------------------------
# PLOTTING
# ------------------------

# Example plot: Fst window counts
plot_Fst <- Fst_window_count_chr %>%
  ggplot(aes(x = Wpos, y = count, color = color_group)) +
  geom_point(size = 2) +
  geom_vline(xintercept = top_window, linetype = "dotted", color = "red", linewidth = 1) +
  scale_x_continuous(breaks = axis_set$center, labels = axis_set$scaf) +
  scale_color_manual(values = color_pairs) +
  labs(
    title = paste0("SexFindR SNP-based analysis for permutation test: ", species),
    subtitle = expression("A. Intersex F"[ST]),
    y = expression("F"[ST] * " outliers per 10k window"),
    x = "Chromosome"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16, angle = 90, vjust = 0.5)
  )

# You can similarly create plots for GWAS, SNP density, nucleotide diversity following this template

# ------------------------
# DISPLAY PLOT
# ------------------------

print(plot_Fst)

# Optionally, save plots to file
# ggsave(filename = file.path(output_dir, paste0(species, "_Fst_plot.png")), plot = plot_Fst, width = 12, height = 8, dpi = 300)

