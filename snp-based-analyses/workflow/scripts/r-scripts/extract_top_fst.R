# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

# Main function
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Check for correct number of arguments
  if (length(args) != 4) {
    stop("Usage: Rscript extract_top_fst_snps.R <fst_file> <ref_file> <species_name> <output_file>")
  }
  
  # Parse arguments
  fst_file <- args[1]
  ref_file <- args[2]
  species_name <- args[3]
  output_file <- args[4]

  # Log processing information
  cat("================================================\n")
  cat("Extracting top FST SNPs for:", species_name, "\n")
  cat("================================================\n")
  cat("Input FST file:       ", fst_file, "\n")
  cat("Reference file:       ", ref_file, "\n")
  cat("Output file:          ", output_file, "\n")
  cat("================================================\n\n")

  # Check if input files exist
  if (!file.exists(fst_file)) {
    stop("FST file not found: ", fst_file)
  }
  
  if (!file.exists(ref_file)) {
    stop("Reference file not found: ", ref_file)
  }

  # Read and preprocess Fst data
  cat("Reading per-SNP FST data...\n")
  Fst_start <- read_tsv(fst_file, show_col_types = FALSE) %>%
    mutate(WEIR_AND_COCKERHAM_FST = ifelse(is.nan(WEIR_AND_COCKERHAM_FST), NA, WEIR_AND_COCKERHAM_FST)) %>%
    replace_na(list(WEIR_AND_COCKERHAM_FST = 0)) %>%
    rename(scaf = CHROM) %>%
    mutate(
      scaf = as.numeric(gsub("Chr", "", scaf)),  # Remove 'Chr' and convert to numeric
      base = as.numeric(POS)
    ) %>% 
    filter(WEIR_AND_COCKERHAM_FST >= 0)
  
  cat("  - Loaded", formatC(nrow(Fst_start), big.mark = ","), "SNPs\n")

  # Load the .fai file to get chromosome lengths
  cat("Reading reference genome index...\n")
  ref <- read.table(ref_file, header = FALSE, 
                    col.names = c("scaf", "LENGTH", "X1", "X2", "X3")) %>%
    mutate(scaf = as.numeric(gsub("Chr", "", scaf))) %>%
    arrange(scaf) %>%
    mutate(cum_start = cumsum(c(0, LENGTH[-n()])))
  
  cat("  - Found", nrow(ref), "chromosomes\n")

  # Merge Fst with reference data
  cat("Merging FST data with reference positions...\n")
  Fst <- left_join(ref, Fst_start, by = "scaf") %>%
    filter(WEIR_AND_COCKERHAM_FST >= 0) %>%
    arrange(scaf, base) %>%
    mutate(Bpos = cum_start + base)

  # Calculate statistics
  cat("\nCalculating FST statistics...\n")
  
  # Top 5% cutoff
  cutoff1 <- round(nrow(Fst) * 0.05)
  Fst_sort_cut1 <- Fst %>%
    arrange(desc(WEIR_AND_COCKERHAM_FST)) %>%
    slice_head(n = cutoff1)
  min_fst_5pct <- min(Fst_sort_cut1$WEIR_AND_COCKERHAM_FST, na.rm = TRUE)
  cat("  - Top 5% threshold:  ", round(min_fst_5pct, 6), 
      " (n=", formatC(cutoff1, big.mark = ","), " SNPs)\n", sep = "")

  # Top 1% cutoff
  cutoff2 <- round(nrow(Fst) * 0.01)
  Fst_sort_cut2 <- Fst %>%
    arrange(desc(WEIR_AND_COCKERHAM_FST)) %>%
    slice_head(n = cutoff2)
  min_fst_1pct <- min(Fst_sort_cut2$WEIR_AND_COCKERHAM_FST, na.rm = TRUE)
  cat("  - Top 1% threshold:  ", round(min_fst_1pct, 6), 
      " (n=", formatC(cutoff2, big.mark = ","), " SNPs)\n", sep = "")

  # Percentile statistics
  quantiles <- quantile(Fst$WEIR_AND_COCKERHAM_FST, c(0.95, 0.99), na.rm = TRUE)
  cat("  - 95th percentile:   ", round(quantiles[1], 6), "\n", sep = "")
  cat("  - 99th percentile:   ", round(quantiles[2], 6), "\n", sep = "")
  cat("  - Mean FST:          ", round(mean(Fst$WEIR_AND_COCKERHAM_FST, na.rm = TRUE), 6), "\n", sep = "")

  # Prepare output: top 5% SNPs with chromosome and position
  cat("\nExtracting top 5% SNPs...\n")
  Fst_sort_cut <- Fst %>%
    arrange(desc(WEIR_AND_COCKERHAM_FST)) %>%
    slice_head(n = cutoff1) %>%
    select(scaf, base)

  # Write output file (no column names, tab-separated)
  write_tsv(Fst_sort_cut, file = output_file, col_names = FALSE)
  cat("  ✓ Top 5% SNPs saved to:", output_file, "\n")
  
  cat("\n================================================\n")
  cat("Top FST SNP extraction completed successfully!\n")
  cat("================================================\n")
}

# Run main
main()