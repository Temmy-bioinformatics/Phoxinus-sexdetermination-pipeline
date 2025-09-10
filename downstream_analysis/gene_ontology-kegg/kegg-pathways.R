# -------------------------
# KEGG Pathway Enrichment
# -------------------------

# Load required libraries
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)      # Change if not human
library(enrichplot)
library(ggplot2)

# -------------------------
# Inputs and directories
# -------------------------
# Use relative paths so it's portable
base_dir <- "Candidate_regions_annotation/RhinePhoxinus/Fst"
input_file <- file.path(base_dir, "Fst_Extracted_Genes_1Mb.tsv")

# -------------------------
# Load and preprocess genes
# -------------------------
gene_data <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE) %>%
  setNames(c("Chromosome", "Start", "End", "Gene_ID", "Gene_Name")) %>%
  distinct() %>%
  filter(Gene_Name != "NA") %>%
  mutate(Gene_Name = toupper(Gene_Name))

# Map gene symbols to Entrez IDs
entrez_ids <- bitr(gene_data$Gene_Name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# -------------------------
# KEGG enrichment analysis
# -------------------------
kegg_enrichment <- enrichKEGG(
  gene = entrez_ids$ENTREZID,
  organism = "hsa",          # Adjust organism code if not human
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Convert results to a data frame with mapped gene symbols
kegg_results_df <- kegg_enrichment@result %>%
  mutate(
    gene_symbols = sapply(geneID, function(ids) {
      genes <- unlist(strsplit(ids, "/"))
      mapped <- entrez_ids %>% filter(ENTREZID %in% genes) %>% pull(SYMBOL)
      paste(mapped, collapse = ", ")
    })
  )

# Filter significant pathways
significant_kegg <- kegg_results_df %>% filter(pvalue < 0.05)

# -------------------------
# Save outputs
# -------------------------
write.table(kegg_results_df, file.path(base_dir, "KEGG_enrichment_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(significant_kegg, file.path(base_dir, "Significant_KEGG_Pathways.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("KEGG pathway analysis complete! Results saved in:", base_dir, "\n")

# -------------------------
# Plotting
# -------------------------
kegg_plot <- ggplot(significant_kegg, aes(x = Count, y = reorder(Description, Count))) +
  geom_point(shape = 21, fill = "black", color = "black", size = 4) +
  geom_text(aes(x = max(Count) + 1, label = gene_symbols), hjust = 0, size = 2.5, color = "black") +
  scale_x_reverse() +
  labs(
    title = expression(paste("KEGG Pathway Enrichment - ", italic("Phoxinus phoxinus"), " SDR regions")),
    x = "Gene Count",
    y = "KEGG Pathways"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = -0.2, face = "bold"),
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12)
  )

# Save plots
ggsave(file.path(base_dir, "KEGG_enrichment_plot_Phoxinus_phoxinus.png"), kegg_plot, width = 11.5, height = 7, dpi = 300)
ggsave(file.path(base_dir, "KEGG_enrichment_plot_Phoxinus_phoxinus.pdf"), kegg_plot, width = 11.5, height = 7)
ggsave(file.path(base_dir, "KEGG_enrichment_plot_Phoxinus_phoxinus.tiff"), kegg_plot, width = 11.5, height = 7, dpi = 300)
ggsave(file.path(base_dir, "KEGG_enrichment_plot_Phoxinus_phoxinus.svg"), kegg_plot, width = 11.5, height = 7)

cat("KEGG enrichment plots saved!\n")
