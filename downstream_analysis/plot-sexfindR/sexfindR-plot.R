# Load required libraries
library(tidyverse)
library(patchwork)
library(ggthemes)

# Clear workspace
rm(list = ls())

# Helper function to calculate cumulative start positions
add_cumulative_start <- function(df) {
  df %>%
    arrange(chr) %>%
    mutate(start = cumsum(Length) - Length)
}

# Set working directory
setwd("C:/Users/TOriowo/Documents/PROJECTS/Phoxinus_SexDetermination_Publication/Haplome_1_results/DanubeCsikii")

###############################################################################
# LOAD DATA

species <- "DanubeCsikii"

# SNP density
SNP <- read_delim(paste0(species, "_SNPdensity_SexFindR.txt")) %>%
  select(scaf, base, mean_MvF_dif, Pvalue) %>%
  mutate(
    scaf = str_remove(scaf, "Chr"),
    base = base + 10000
  ) %>%
  filter(Pvalue <= 0.05)

# Rank SNP density windows
SNP_filter_arrange <- SNP %>%
  arrange(Pvalue, desc(abs(mean_MvF_dif))) %>%
  rownames_to_column(var = "SNPdensity_rank") %>%
  mutate(
    scaf = as.numeric(scaf),
    SNPdensity_rank = as.numeric(SNPdensity_rank)
  )

# GWAS window counts
GWAS_window_count <- read_delim(paste0(species, "_gwas_window_count.txt")) %>%
  arrange(desc(count)) %>%
  rownames_to_column(var = "GWAScount_rank") %>%
  mutate(
    scaf = as.numeric(scaf),
    GWAScount_rank = as.numeric(GWAScount_rank)
  )

# Fst window counts
Fst_window_count <- read_delim(paste0(species, "_fst_window_count.txt")) %>%
  arrange(desc(count)) %>%
  rownames_to_column(var = "Fstcount_rank") %>%
  mutate(Fstcount_rank = as.numeric(Fstcount_rank))

# Pi differences
pi <- read_delim(paste0(species, "_pi_10kb_sex_differences.txt")) %>%
  select(scaf, BIN_START, pi_dif_M_minus_F) %>%
  arrange(desc(abs(pi_dif_M_minus_F))) %>%
  rownames_to_column(var = "pi_rank") %>%
  rename(base = BIN_START) %>%
  mutate(
    pi_rank = as.numeric(pi_rank),
    base = base + 9999,
    scaf = as.numeric(scaf)
  )

# Reference genome scaffold lengths and cumulative start positions
ref_file <- "../Final_Phoxy_1_rearranged.fasta.fai"
ref <- read.table(ref_file, header = FALSE, col.names = c("scaf", "LENGTH", "X1", "X2", "X3")) %>%
  mutate(scaf = as.numeric(str_remove(scaf, "Chr"))) %>%
  arrange(scaf) %>%
  mutate(cum_start = cumsum(c(0, head(LENGTH, -1))))

# Assign alternating color groups by chromosome
color_pairs <- c("#FC766AFF", "#5B84B1FF")
ref <- ref %>%
  mutate(color_group = color_pairs[(scaf %% 2) + 1])

# Calculate chromosome midpoints for x-axis labels
axis.set <- ref %>%
  group_by(scaf) %>%
  summarize(center = (cum_start + LENGTH / 2), .groups = "drop")

# Join datasets with reference info
Fst_window_count_chr <- left_join(Fst_window_count, ref, by = "scaf")
GWAS_window_count_chr <- left_join(GWAS_window_count, ref, by = "scaf")
SNP_filter_arrange_chr <- left_join(SNP_filter_arrange, ref, by = "scaf")

# Combine all SNP, Fst, GWAS data
SNP_Fst <- full_join(SNP_filter_arrange_chr, Fst_window_count_chr, by = c("scaf", "base"))
SNP_Fst_GWAS <- left_join(SNP_Fst, GWAS_window_count, by = c("scaf", "base"))

all <- full_join(SNP_Fst_GWAS, pi, by = c("scaf", "base")) %>%
  select(-count.x, -count.y) %>%
  mutate(total_rank = GWAScount_rank + SNPdensity_rank + Fstcount_rank + pi_rank)

###############################################################################
# SNP counts per window and average
snp_counts_final <- all %>%
  group_by(scaf, base) %>%
  summarise(snp_count = n(), .groups = "drop")

average_snps_per_window_final <- mean(snp_counts_final$snp_count)

write_tsv(snp_counts_final, paste0(species, "_Final_SNPs_per_window.txt"))
write_lines(paste("Average SNPs per window:", average_snps_per_window_final),
            paste0(species, "_Final_Average_SNPs_per_window.txt"))

###############################################################################
# Spearman correlations between rankings
correlations <- c(
  "SNP_density vs GWAS_count" = cor.test(all$SNPdensity_rank, all$GWAScount_rank, method = 'spearman')$estimate,
  "SNP_density vs Fst_count" = cor.test(all$SNPdensity_rank, all$Fstcount_rank, method = 'spearman')$estimate,
  "SNP_density vs pi" = cor.test(all$SNPdensity_rank, all$pi_rank, method = 'spearman')$estimate,
  "GWAS_count vs Fst_count" = cor.test(all$GWAScount_rank, all$Fstcount_rank, method = 'spearman')$estimate,
  "GWAS_count vs pi" = cor.test(all$GWAScount_rank, all$pi_rank, method = 'spearman')$estimate,
  "Fst_count vs pi" = cor.test(all$Fstcount_rank, all$pi_rank, method = 'spearman')$estimate
)

# Sort by absolute correlation strength descending
sorted_correlations <- sort(abs(correlations), decreasing = TRUE)
names(sorted_correlations) <- names(correlations)[order(abs(correlations), decreasing = TRUE)]

# Write correlations to file
sink(paste0(species, "_correlation_rankings.txt"))
cat("Ranked Spearman's Correlations (from strongest to weakest):\n\n")
for (i in seq_along(sorted_correlations)) {
  cat(i, ": ", names(sorted_correlations)[i], " - rho = ", correlations[names(sorted_correlations)[i]], "\n", sep = "")
}
sink()

###############################################################################
# Filter top-ranked windows
top_all <- all %>%
  filter(
    GWAScount_rank <= 100,
    Fstcount_rank <= 100,
    pi_rank <= 100,
    SNPdensity_rank <= 100
  )

# 3-way combinations
top_rank_SFP <- all %>% filter(SNPdensity_rank <= 100, Fstcount_rank <= 100, pi_rank <= 100)
top_rank_GFP <- all %>% filter(GWAScount_rank <= 100, Fstcount_rank <= 100, pi_rank <= 100)
top_rank_GFS <- all %>% filter(GWAScount_rank <= 100, Fstcount_rank <= 100, SNPdensity_rank <= 100)
top_rank_SGP <- all %>% filter(SNPdensity_rank <= 100, GWAScount_rank <= 100, pi_rank <= 100)

# 2-way combinations
top_rank_SP <- all %>% filter(SNPdensity_rank <= 100, pi_rank <= 100)
top_rank_GF <- all %>% filter(GWAScount_rank <= 100, Fstcount_rank <= 100)
top_rank_GS <- all %>% filter(GWAScount_rank <= 100, SNPdensity_rank <= 100)
top_rank_SF <- all %>% filter(SNPdensity_rank <= 100, Fstcount_rank <= 100)
top_rank_GP <- all %>% filter(GWAScount_rank <= 100, pi_rank <= 100)
top_rank_FP <- all %>% filter(Fstcount_rank <= 100, pi_rank <= 100)

# Save filtered datasets (dropping total_rank for clarity)
write_tsv(select(top_all, -total_rank), "SexFindR_candidate_windows_SNP_GWAS_Fst_pi.txt")
write_tsv(select(top_rank_GFP, -total_rank), "SexFindR_candidate_windows_GWAS_Pi_Fst.txt")
write_tsv(select(top_rank_SFP, -total_rank), "SexFindR_candidate_windows_SNP_Fst_Pi.txt")
write_tsv(select(top_rank_GFS, -total_rank), "SexFindR_candidate_windows_GWAS_Fst_SNP.txt")
write_tsv(select(top_rank_SGP, -total_rank), "SexFindR_candidate_windows_SNP_GWAS_Pi.txt")

write_tsv(select(top_rank_SP, -total_rank), "SexFindR_candidate_windows_SNP_Pi.txt")
write_tsv(select(top_rank_GF, -total_rank), "SexFindR_candidate_windows_GWAS_Fst.txt")
write_tsv(select(top_rank_GS, -total_rank), "SexFindR_candidate_windows_GWAS_SNP.txt")
write_tsv(select(top_rank_SF, -total_rank), "SexFindR_candidate_windows_SNP_Fst.txt")
write_tsv(select(top_rank_GP, -total_rank), "SexFindR_candidate_windows_GWAS_Pi.txt")
write_tsv(select(top_rank_FP, -total_rank), "SexFindR_candidate_windows_Fst_Pi.txt")

###############################################################################
# Calculate window positions for highlighting in plots
calc_Wpos <- function(df) {
  df %>% mutate(Wpos = cum_start + base)
}

top_points_list <- list(
  top_all, top_rank_GFP, top_rank_SFP, top_rank_GFS, top_rank_SGP,
  top_rank_SP, top_rank_SF, top_rank_GF, top_rank_GS, top_rank_GP, top_rank_FP
)

names(top_points_list) <- c("top_points", "top_points_GFP", "top_points_SFP", "top_points_GFS", "top_points_SGP",
                            "top_points_SP", "top_points_SF", "top_points_GF", "top_points_GS", "top_points_GP", "top_points_FP")

for (nm in names(top_points_list)) {
  assign(nm, calc_Wpos(top_points_list[[nm]]))
  assign(sub("points", "window", nm), get(nm)$Wpos)
}

###############################################################################
# Plotting

fullname <- "Phoxinus phoxinus"

# Fst plot
Fst_window_count_chr <- Fst_window_count_chr %>% mutate(Wpos = cum_start + base)
a <- ggplot(Fst_window_count_chr, aes(x = Wpos, y = count)) +
  geom_point(aes(color = color_group), size = 2) +
  geom_vline(xintercept = top_window, linetype = "dotted", color = "red", linewidth = 1) +
  scale_x_continuous(breaks = axis.set$center, labels = axis.set$scaf) +
  scale_color_manual(values = color_pairs, guide = "none") +
  scale_y_continuous(limits = c(0, max(Fst_window_count_chr$count) * 1.05), expand = c(0, 0)) +
  labs(
    title = bquote("SexFindR SNP-based analysis for permutation test: "~italic(.(fullname))),
    subtitle = expression("A. Intersex F"[ST]),
    y = expression("F"[ST]*" outliers per 10k window"),
    x = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 22, angle = 90, hjust = 0.5, vjust = 0.5),
    plot.title = element_text(size = 30)
  )

# GWAS plot
GWAS_window_count_chr <- GWAS_window_count_chr %>% mutate(Wpos = cum_start + base)
b <- ggplot(GWAS_window_count_chr, aes(x = Wpos, y = count)) +
  geom_point(aes(color = color_group)) +
  geom_vline(xintercept = top_window, linetype = "dotted", color = "red", linewidth = 1) +
  scale_x_continuous(breaks = axis.set$center, labels = axis.set$scaf) +
  scale_color_manual(values = color_pairs, guide = "none") +
  scale_y_continuous(limits = c(0, max(GWAS_window_count_chr$count) * 1.05), expand = c(0, 0)) +
  labs(
    title = "B. GWAS of sex",
    y = "GWAS outliers per 10k window",
    x = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 25, angle = 90, hjust = 0.5, vjust = 0.5),
    plot.title = element_text(size = 30)
  )

# SNP density plot
SNP_filter_arrange_chr <- SNP_filter_arrange_chr %>% mutate(Wpos = cum_start + base)
c <- ggplot(SNP_filter_arrange_chr, aes(x = Wpos, y = mean_MvF_dif)) +
  geom_point(aes(color = color_group)) +
  geom_vline(xintercept = top_window, linetype = "dotted", color = "red", linewidth = 1) +
  scale_x_continuous(breaks = axis.set$center, labels = axis.set$scaf) +
  scale_color_manual(values = color_pairs, guide = "none") +
  scale_y_continuous(
    limits = c(min(SNP_filter_arrange_chr$mean_MvF_dif, na.rm = TRUE) * 1.25,
               max(SNP_filter_arrange_chr$mean_MvF_dif, na.rm = TRUE) * 1.25),
    expand = c(0, 0)
  ) +
  labs(
    title = "C. SNP density",
    y = "10kb SNP Density \n(male mean - female mean)",
    x = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 25, angle = 90, hjust = 0.5, vjust = 0.5),
    plot.title = element_text(size = 30)
  )

# Nucleotide diversity plot
pi_chr <- left_join(pi, ref, by = "scaf") %>% mutate(Wpos = cum_start + base)
d <- ggplot(pi_chr, aes(x = Wpos, y = pi_dif_M_minus_F)) +
  geom_point(aes(color = color_group)) +
  geom_vline(xintercept = top_window, linetype = "dotted", color = "red", linewidth = 1) +
  scale_x_continuous(breaks = axis.set$center, labels = axis.set$scaf) +
  scale_color_manual(values = color_pairs, guide = "none") +
  scale_y_continuous(
    limits = c(min(pi_chr$pi_dif_M_minus_F, na.rm = TRUE) * 1.25,
               max(pi_chr$pi_dif_M_minus_F, na.rm = TRUE) * 1.25),
    expand = c(0, 0)
  ) +
  labs(
    title = "D. Nucleotide diversity",
    y = "10kb \u03C0 \n(male mean - female mean)",
    x = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 20, vjust = 0.5),
    axis.title.x = element_text(size = 25),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 25, angle = 90, hjust = 0.5, vjust = 0.5),
    plot.title = element_text(size = 30)
  )

# Combine plots
combined_plot <- a / b / c / d
combined_plot

# Uncomment to save plots
# ggsave(paste0("Results_plots/", species, "_SexFindR_top_combined_results.pdf"), combined_plot, width = 23, height = 27)
# ggsave(paste0("Results_plots/", species, "_SexFindR_top_combined_results.tiff"), combined_plot, width = 23, height = 27, dpi = 300)
# ggsave(paste0("Results_plots/", species, "_SexFindR_combined_top_results.png"), combined_plot, width = 23, height = 27, dpi = 600)


#Display the final plot
#print(combined_plot)


#################### PER CHROMOSOME PLOT######################################################## # nolint
###############################################################################################
# # Dynamically set the linkage group
# lg <- "3"  # You can change this to any desired linkage group
# 
# # Extract the chromosome end value from the LENGTH column
# chromosome_end <- top_all %>%
#   filter(scaf == lg) %>%
#   pull(LENGTH) %>%
#   unique()
# 
# x_scale_limits <- c(0, chromosome_end)  # X-axis is now between 0 and the LENGTH value
# 
# # Ensure your top_window_lg stays consistent
# top_points_lg <- top_all %>%
#   filter(scaf == lg) %>%
#   mutate(Wpos = base)
# 
# top_window_lg <- top_points_lg$Wpos
# 
# # Modify the first plot (e) to include the gene track
# e <- Fst_window_count_chr %>%
#   filter(scaf == lg) %>%
#   ggplot(aes(x = base, y = count)) +
#   geom_point(aes(fill = color_group), color = "#66A61E", size = 1) +
#   labs(
#     x = NULL,
#     y = expression("F"[ST] * " outliers per 10k window"),
#     title = bquote("SexFindR combined SNP-based analysis for" ~ italic(.(fullname)) ~ " chromosome " ~ .(lg)),
#     subtitle = expression("A. Intersex F"[ST])
#   ) +
#   gene_track +  # Add gene track
#   geom_vline(xintercept = top_window_lg, linetype = "dotted", color = "red", linewidth = 1) +
#   scale_x_continuous(limits = x_scale_limits, expand = c(0, 0), labels = scales::comma) +
#   theme_bw() +
#   theme(
#     legend.position = "none",
#     axis.text.x = element_blank(),
#     axis.title.x = element_blank(),
#     axis.text.y = element_text(size = 11),
#     axis.title.y = element_text(size = 15, angle = 90, hjust = 0.5, vjust = 0.5, face = "plain"),
#     title = element_text(size = 15),
#     plot.title = element_text(size = 20)
#   )
# 
# # The other plots (f, g, h) without the gene track
# f <- GWAS_window_count_chr %>%
#   filter(scaf == lg) %>%
#   ggplot(aes(x = base, y = count)) +
#   geom_point(aes(fill = color_group), color = "#66A61E", size = 1) +
#   labs(
#     x = NULL,
#     y = "GWAS outliers per 10kb window",
#     subtitle = "B. GWAS"
#   ) +
#   geom_vline(xintercept = top_window_lg, linetype = "dotted", color = "red", linewidth = 1) +
#   scale_x_continuous(limits = x_scale_limits, expand = c(0, 0), labels = scales::comma) +
#   theme_bw() +
#   theme(
#     legend.position = "none",
#     axis.text.x = element_blank(),
#     axis.title.x = element_blank(),
#     axis.text.y = element_text(size = 11),
#     axis.title.y = element_text(size = 15, angle = 90, hjust = 0.5, vjust = 0.5, face = "plain"),
#     title = element_text(size = 15),
#     plot.title = element_text(size = 20)
#   )
# 
# g <- SNP_filter_arrange_chr %>%
#   filter(scaf == lg) %>%
#   ggplot(aes(x = base, y = mean_MvF_dif)) +
#   geom_point(aes(fill = color_group), color = "#66A61E", size = 1) +
#   geom_vline(xintercept = top_window_lg, linetype = "dotted", color = "red", linewidth = 1) +
#   labs(
#     x = NULL,
#     y = "10kb SNP density \n(male mean - female mean)",
#     subtitle = "C. SNP density"
#   ) +
#   scale_x_continuous(limits = x_scale_limits, expand = c(0, 0), labels = scales::comma) +
#   theme_bw() +
#   theme(
#     legend.position = "none",
#     axis.text.x = element_blank(),
#     axis.title.x = element_blank(),
#     axis.text.y = element_text(size = 11),
#     axis.title.y = element_text(size = 15, angle = 90, hjust = 0.5, vjust = 0.5, face = "plain"),
#     title = element_text(size = 15),
#     plot.title = element_text(size = 20)
#   )
# 
# h <- pi_chr %>%
#   filter(scaf == lg) %>%
#   ggplot(aes(x = base, y = pi_dif_M_minus_F)) +
#   geom_point(aes(fill = color_group), color = "#66A61E", size = 1) +
#   geom_vline(xintercept = top_window_lg, linetype = "dotted", color = "red", linewidth = 1) +
#   labs(
#     x = "Genome position (bp)",
#     y = "10kb \u03C0 \n(male mean - female mean)",
#     subtitle = "D. Nucleotide Diversity"
#   ) +
#   scale_x_continuous(limits = x_scale_limits, expand = c(0, 0), labels = scales::comma) +
#   theme_bw() +
#   theme(
#     legend.position = "none",
#     axis.text.x = element_text(size = 12),
#     axis.title.x = element_text(size = 15),
#     axis.text.y = element_text(size = 11),
#     axis.title.y = element_text(size = 15, angle = 90, hjust = 0.5, vjust = 0.5, face = "plain"),
#     title = element_text(size = 15),
#     plot.title = element_text(size = 20)
#   )
# 
# # Combine all plots
# combined_plot <- (e / f / g / h)
# 
# # Save the plot as PDF, TIFF, and PNG
# output_dir <- "Results_plots"
# dir.create(output_dir, showWarnings = FALSE)
# pdf_filename <- file.path(output_dir, paste0(fullname, "_LG", lg, "_combined_analysis.pdf"))
# tiff_filename <- file.path(output_dir, paste0(fullname, "_LG", lg, "_combined_analysis.tiff"))
# png_filename <- file.path(output_dir, paste0(fullname, "_LG", lg, "_combined_analysis.png"))
# 
# # Save the combined plot
# ggsave(pdf_filename, plot = combined_plot, width = 14, height = 15, units = "in")
# ggsave(tiff_filename, plot = combined_plot, width = 14, height = 15, units = "in", dpi = 300, device = "tiff")
# ggsave(png_filename, plot = combined_plot, width = 14, height = 15, units = "in", dpi = 600, device = "png")
