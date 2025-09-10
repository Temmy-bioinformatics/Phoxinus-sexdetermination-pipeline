# -------------------------
# GO Enrichment for Flanking Genes
# -------------------------

# Load required libraries
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)      # Change if not using human annotation
library(enrichplot)

# -------------------------
# Settings
# -------------------------
species <- "DanubeCsikii"
species_name <- "Phoxinus csikii"
working_dir <- file.path("Candidate_regions_annotation", species, "select_regions")
setwd(working_dir)

# Get all TSV files that match the flanking gene pattern
tsv_files <- list.files(pattern = "^flanking_genes_\\d+_\\d+_\\d+_\\d+\\.tsv$")

# -------------------------
# Loop through each file
# -------------------------
for (file_path in tsv_files) {
  cat("Processing file:", file_path, "\n")

  region_id <- tools::file_path_sans_ext(basename(file_path))

  # Load gene data
  gene_data <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE) %>%
    setNames(c("Chromosome", "Start", "End", "Gene_ID", "Gene_Name", "SV")) %>%
    distinct() %>%
    filter(Gene_Name != "NA") %>%
    mutate(Gene_Name = toupper(Gene_Name))

  # Save gene names for external use
  write.table(
    gene_data$Gene_Name,
    paste0(region_id, "_500kb_flanking_gene_list.txt"),
    quote = FALSE, row.names = FALSE, col.names = FALSE
  )

  # Convert gene symbols to Entrez IDs
  entrez_ids <- bitr(gene_data$Gene_Name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

  # Run GO enrichment (Biological Processes)
  go_enrichment <- enrichGO(
    gene = gene_data$Gene_Name,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )

  # Save complete GO enrichment results
  go_results_df <- go_enrichment@result
  write.table(go_results_df, paste0(region_id, "_500kb_GO_enrichment_full.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  # Filter significant GO terms (p < 0.01)
  significant_go <- go_results_df %>% filter(pvalue < 0.01)

  if (nrow(significant_go) > 0) {
    # Map genes to GO terms
    significant_go$gene_symbols <- sapply(significant_go$ID, function(go_id) {
      genes <- go_enrichment@result[go_enrichment@result$ID == go_id, "geneID"]
      paste(unlist(strsplit(genes, "/")), collapse = ", ")
    })

    # Save significant terms
    write.table(significant_go, paste0(region_id, "_500kb_GO_terms_p0.01.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    # Create black & white GO bubble plot
    go_plot <- ggplot(significant_go, aes(x = Count, y = reorder(Description, Count))) +
      geom_point(shape = 21, fill = "black", color = "black", size = 4) +
      geom_text(aes(x = max(Count) + 1, label = gene_symbols), hjust = 0, size = 2.5, color = "black") +
      scale_x_reverse(breaks = seq(0, floor(max(significant_go$Count)), by = 1)) +
      labs(
        title = paste("GO: Biological Processes ±500 kb from SDR in", species_name, "\n(", region_id, ")"),
        x = "Gene Count",
        y = "GO Terms"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 14, hjust = 0),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13)
      )

    # Export plot in multiple formats
    for (ext in c("png", "pdf", "tiff", "svg")) {
      ggsave(
        filename = paste0(region_id, "_GO_BP_500kb_plot.", ext),
        plot = go_plot,
        width = 12,
        height = 10,
        dpi = 600,
        units = "in"
      )
    }
  } else {
    cat("No significant GO terms (p < 0.01) for:", region_id, "\n")
  }

  cat("Finished processing", file_path, "\n\n")
}
