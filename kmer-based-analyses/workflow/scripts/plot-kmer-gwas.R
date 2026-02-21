library(tidyverse)
library(cowplot)
library(scales)

# Main function
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Check for correct number of arguments
  if (length(args) != 5) {
    stop(
      "Usage: Rscript plot_kmer_blast.R ",
      "<blast_file> <ref_file> <species_name> <output_dir> <species_code>"
    )
  }
  
  # Parse arguments
  blast_file <- args[1]
  ref_file <- args[2]
  species_name <- args[3]  # e.g., "Phoxinus phoxinus"
  output_dir <- args[4]
  species_code <- args[5]  # e.g., "RhinePhoxinus"
  
  # Log processing information
  cat("================================================\n")
  cat("K-mer BLAST Plotting for:", species_name, "\n")
  cat("================================================\n")
  cat("BLAST file:           ", blast_file, "\n")
  cat("Reference file:       ", ref_file, "\n")
  cat("Output directory:     ", output_dir, "\n")
  cat("Species code:         ", species_code, "\n")
  cat("================================================\n\n")
  
  # Check if input files exist
  if (!file.exists(blast_file)) stop("BLAST file not found: ", blast_file)
  if (!file.exists(ref_file)) stop("Reference file not found: ", ref_file)
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Define color pairs for species
  color_pairs <- list(
    "DanubeCsikii" = c("#9897A9", "#D95E02"),
    "RhinePhoxinus" = c("#9897A9", "#66A61E"),
    "SeinePhoxinus" = c("#9897A9", "#ab9519")
  )
  
  # P-values for species
  p_values <- list(
    "RhinePhoxinus" = "1.183 × 10⁻¹¹",
    "DanubeCsikii" = "1.542 × 10⁻⁰⁸",
    "SeinePhoxinus" = "TBD"
  )
  
  # Get colors
  cols <- if (species_code %in% names(color_pairs)) {
    color_pairs[[species_code]]
  } else {
    c("#9897A9", "#D95E02")
  }
  
  # Get p-value
  p_value <- if (species_code %in% names(p_values)) {
    p_values[[species_code]]
  } else {
    "Not available"
  }
  
  cat("Using colors:", cols[1], ",", cols[2], "\n")
  cat("P-value:", p_value, "\n\n")
  
  # ------------------------
  # LOAD AND PROCESS DATA
  # ------------------------
  
  # Load reference file
  cat("Loading reference genome index...\n")
  ref <- read.table(ref_file, header = FALSE,
                    col.names = c("scaf", "LENGTH", "X1", "X2", "X3")) %>%
    mutate(scaf = as.numeric(gsub("Chr", "", scaf))) %>%
    arrange(scaf)
  
  cat("  - Found", nrow(ref), "chromosomes\n")
  
  # Load BLAST file
  cat("Loading BLAST results...\n")
  blast <- read_tsv(
    blast_file,
    col_names = c(
      "query_sid", "scaf", "ident.match_p", "alignmt_len", "mismatch_n",
      "gapopen_n", "query_start", "query_end", "ref_start", "ref_end",
      "exp_value", "bitscore"
    ),
    col_types = cols(
      query_sid = col_character(),
      scaf = col_character(),
      ident.match_p = col_double(),
      exp_value = col_double()
    ),
    show_col_types = FALSE
  ) %>%
    mutate(scaf = as.numeric(str_remove(scaf, "^Chr"))) %>%
    filter(!is.na(exp_value))
  
  cat("  - Loaded", nrow(blast), "BLAST hits\n")
  
  # Best alignment per k-mer
  cat("Filtering best alignments per k-mer...\n")
  blast_filter <- blast %>%
    group_by(query_sid) %>%
    slice_min(order_by = exp_value, n = 1) %>%
    ungroup()
  
  cat("  - Filtered to", nrow(blast_filter), "k-mers\n")
  
  # Save summary
  cat("Writing mapping summary...\n")
  summary_row <- tibble(
    species = species_code,
    total_kmers_searched = length(unique(blast$query_sid)),
    kmers_with_blast_hits = nrow(blast_filter)
  )
  
  summary_file <- file.path(
    output_dir,
    paste0(species_code, "_kmer_blast_mapping_summary.csv")
  )
  write_csv(summary_row, summary_file)
  cat("  ✓ Summary saved to:", summary_file, "\n")
  
  # Count hits per chromosome
  cat("Counting hits per chromosome...\n")
  blast_chr <- blast_filter %>%
    count(scaf) %>%
    complete(scaf = ref$scaf, fill = list(n = 0)) %>%
    left_join(ref, by = "scaf") %>%
    arrange(scaf) %>%
    mutate(prop_n = n / LENGTH)
  
  # ------------------------
  # GENERATE PLOT
  # ------------------------
  
  cat("Generating chromosome-wide plot...\n")
  
  chr_plot <- ggplot(blast_chr,
                     aes(x = scaf, y = prop_n * 1e6,
                         fill = factor(scaf %% 2))) +
    geom_col() +
    theme_bw() +
    labs(
      title = bquote(
        "Distribution of top k-mer contigs with BLAST hits in"
        ~ italic(.(species_name))
      ),
      subtitle = paste(
        "Mapped k-mers per chromosome, normalized per Mb (lowest p-value =",
        p_value, ")"
      ),
      x = "Chromosome",
      y = "K-mer derived contigs per Mb"
    ) +
    scale_x_continuous(breaks = ref$scaf, labels = ref$scaf) +
    scale_fill_manual(values = cols) +
    theme(
      axis.text.x = element_text(size = 12, vjust = 0.5),
      axis.title.x = element_text(size = 14),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.position = "none"
    )
  
  # ------------------------
  # SAVE PLOT
  # ------------------------
  
  cat("Saving plots...\n")
  
  ggsave(file.path(output_dir,
                   paste0(species_code, "_kmer_counts_per_chr.png")),
         plot = chr_plot, width = 11.5, height = 7, dpi = 300)
  
  ggsave(file.path(output_dir,
                   paste0(species_code, "_kmer_counts_per_chr.pdf")),
         plot = chr_plot, width = 11.5, height = 7)
  
  ggsave(file.path(output_dir,
                   paste0(species_code, "_kmer_counts_per_chr.svg")),
         plot = chr_plot, width = 11.5, height = 7)
  
  cat("  ✓ Plots saved (PNG, PDF, SVG)\n")
  
  cat("\n================================================\n")
  cat("K-mer BLAST plotting completed successfully!\n")
  cat("================================================\n")
}

# Run main
main()
