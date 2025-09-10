# Load required libraries
library(tidyverse)
library(cowplot)
library(scales)

# Set project working directory
setwd("kmerGWAS")

# Define color pairs for species plots
color_pairs <- list(
  "DanubeCsikii" = c("#9897A9", "#D95E02"),
  "RhinePhoxinus" = c("#9897A9", "#66A61E")
)

# Function to generate species plots and write summary counts
generate_species_plots <- function(species, species_name, color_pairs, ref_file, normal_blast_file, top_blast_file, chr_zoom = 3) {
  # Load reference genome index (.fai file)
  ref <- read.table(ref_file, header = FALSE, col.names = c("scaf", "LENGTH", "X1", "X2", "X3")) %>%
    mutate(scaf = as.numeric(gsub("Chr", "", scaf))) %>%
    arrange(scaf)
  
  # Helper to load and process BLAST hits
  load_blast <- function(file) {
    read_tsv(
      file,
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
      )
    ) %>%
      mutate(scaf = as.numeric(str_remove(scaf, "^Chr"))) %>%
      filter(!is.na(exp_value))
  }
  
  normal_blast <- load_blast(normal_blast_file)
  top_blast <- load_blast(top_blast_file)
  
  # Keep the best alignment per k-mer (lowest e-value)
  normal_filter <- normal_blast %>%
    group_by(query_sid) %>%
    slice_min(order_by = exp_value, n = 1) %>%
    ungroup()
  top_filter <- top_blast %>%
    group_by(query_sid) %>%
    slice_min(order_by = exp_value, n = 1) %>%
    ungroup()
  
  # Save summary table
  summary_row <- tibble(
    species = species,
    total_kmers_searched = length(unique(normal_blast$query_sid)),
    kmers_with_blast_hits = nrow(normal_filter),
    top001_kmers_searched = length(unique(top_blast$query_sid)),
    top001_kmers_with_hits = nrow(top_filter)
  )
  write_csv(summary_row, file = "Results_plots/kmer_blast_mapping_summary.csv", append = file.exists("Results_plots/kmer_blast_mapping_summary.csv"))
  
  # Compute genomic positions for mapped k-mers
  normal_pos <- normal_filter %>%
    left_join(ref, by = "scaf") %>%
    arrange(scaf, ref_start) %>%
    mutate(Kpos = ((query_end - query_start) / 2) + ref_start,
           fill = ifelse(scaf %% 2 == 0, "Even", "Odd"))
  
  top_pos <- top_filter %>%
    left_join(ref, by = "scaf") %>%
    arrange(scaf, ref_start) %>%
    mutate(Kpos = ((query_end - query_start) / 2) + ref_start,
           fill = ifelse(scaf %% 2 == 0, "Even", "Odd"))
  
  # Count hits per chromosome and normalize by chromosome length
  normal_chr <- normal_filter %>%
    count(scaf) %>%
    complete(scaf = ref$scaf, fill = list(n = 0)) %>%
    left_join(ref, by = "scaf") %>%
    arrange(scaf) %>%
    mutate(prop_n = n / LENGTH)
  
  top_chr <- top_filter %>%
    count(scaf) %>%
    complete(scaf = ref$scaf, fill = list(n = 0)) %>%
    left_join(ref, by = "scaf") %>%
    arrange(scaf) %>%
    mutate(prop_n = n / LENGTH)
  
  # Colors
  cols <- color_pairs[[species]]
  alt_color <- cols[2]
  
  # Frame A: all k-mer hits per chromosome
  chr_plot <- ggplot(normal_chr, aes(x = scaf, y = prop_n * 1e6, fill = factor(scaf %% 2))) +
    geom_col() +
    theme_bw() +
    labs(
      title = bquote("Distribution of k-mer BLAST hits in" ~ italic(.(species_name))),
      x = "Chromosome",
      y = "Hits per Mb"
    ) +
    scale_x_continuous(breaks = ref$scaf, labels = ref$scaf) +
    scale_fill_manual(values = cols) +
    theme(legend.position = "none")
  
  # Frame C: top 0.001% k-mer hits per chromosome
  top_chr_plot <- ggplot(top_chr, aes(x = scaf, y = prop_n * 1e6, fill = factor(scaf %% 2))) +
    geom_col() +
    theme_bw() +
    labs(
      title = bquote("Top 0.001% k-mer BLAST hits in" ~ italic(.(species_name))),
      x = "Chromosome",
      y = "Top 0.001% hits per Mb"
    ) +
    scale_x_continuous(breaks = ref$scaf, labels = ref$scaf) +
    scale_fill_manual(values = cols) +
    theme(legend.position = "none")
  
  # Frame B: zoom on chr3, top 0.001% hits
  zoom_top_chr3_plot <- top_pos %>%
    filter(scaf == chr_zoom) %>%
    ggplot(aes(x = Kpos, y = ident.match_p)) +
    geom_point(size = 2, alpha = 0.7, color = alt_color) +
    scale_x_continuous(labels = comma) +
    labs(
      title = bquote("Top 0.001% k-mer hits on Chr" ~ .(chr_zoom) ~ "in" ~ italic(.(species_name))),
      x = "Genomic position",
      y = "BLAST similarity (%)"
    ) +
    theme_bw()
  
  # Frame D: zoom on chr3, all k-mers
  zoom_plot <- normal_pos %>%
    filter(scaf == chr_zoom) %>%
    ggplot(aes(x = Kpos, y = ident.match_p, color = fill)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_manual(values = "#D95E02", guide = "none") +
    scale_x_continuous(labels = comma) +
    labs(
      title = bquote("All k-mer hits on Chr" ~ .(chr_zoom) ~ "in" ~ italic(.(species_name))),
      x = "Genomic position",
      y = "BLAST similarity (%)"
    ) +
    theme_bw()
  
  # Save plots to results directory
  ggsave(paste0("Results_plots/", species, "_frame_A_all_hits.png"), chr_plot, width = 11.5, height = 7)
  ggsave(paste0("Results_plots/", species, "_frame_B_zoom_chr", chr_zoom, "_top001.png"), zoom_top_chr3_plot, width = 11.5, height = 7)
  ggsave(paste0("Results_plots/", species, "_frame_C_top001_hits.png"), top_chr_plot, width = 11.5, height = 7)
  ggsave(paste0("Results_plots/", species, "_frame_D_zoom_chr", chr_zoom, "_all_hits.png"), zoom_plot, width = 11.5, height = 7)
  
  list(chr_plot = chr_plot, top_chr_plot = top_chr_plot,
       zoom_top_chr3_plot = zoom_top_chr3_plot, zoom_plot = zoom_plot)
}

# Define species datasets
species_data <- list(
  list(species = "RhinePhoxinus", species_name = "Phoxinus phoxinus"),
  list(species = "DanubeCsikii", species_name = "Phoxinus csikii")
)

# Reference genome index file (relative path)
ref_file <- "reference/Final_Phoxy_1_rearranged.fasta.fai"

# Generate and store plots for each species
plots <- lapply(species_data, function(data) {
  generate_species_plots(
    species = data$species,
    species_name = data$species_name,
    color_pairs = color_pairs,
    ref_file = ref_file,
    normal_blast_file = paste0(data$species, "/kmers_blast.out"),
    top_blast_file = paste0(data$species, "/top0.001/kmers_blast.out")
  )
})

# Combine and save final multi-panel figure
final_combined_plot <- plot_grid(
  plots[[1]]$chr_plot,
  plots[[1]]$zoom_top_chr3_plot,
  plots[[2]]$chr_plot,
  plots[[2]]$zoom_plot,
  labels = c("A.", "B.", "C.", "D."),
  ncol = 2
)
ggsave("Results_plots/final_combined_kmer_blast_plot.png", final_combined_plot, width = 24, height = 14)
