# Load necessary libraries
library(gplots)
library(pheatmap)
library(vcfR)
library(reshape2)

# Define species (adjust for each run)
Species <- "RhinePhoxinus"
Species_name <- "Phoxinus phoxinus"

# Set working directory relative to project root
base_dir <- file.path("heterozygote_plots", Species)
setwd(base_dir)

# Load phenotype info (no header assumed)
pheno_file <- paste0(Species, ".pheno.txt")
FullSpeciesInfo <- read.delim(pheno_file, header = FALSE, stringsAsFactors = FALSE)
colnames(FullSpeciesInfo) <- c("IDs", "Species", "Sex", "Drainage", "SampleName")

# Convert sex codes to "Male" and "Female"
FullSpeciesInfo$Sex <- ifelse(FullSpeciesInfo$Sex %in% c("1", 1, "M", "Male"), "Male",
                              ifelse(FullSpeciesInfo$Sex %in% c("2", 2, "F", "Female"), "Female", NA))

# Heatmap colors and sex annotations
my_palette <- colorRampPalette(c("#B0B0B0", "#0E3F5C", "#2A6D7A", "#8FCCB4"))(4)
annotation_colors <- list(Sex = c(Male = "#4D83C5", Female = "#D62221"))

# Create output directories for plots
dirs <- c("pdf", "png", "tiff")
for (d in dirs) if (!dir.exists(d)) dir.create(d)

# Find all VCF files (.vcf.gz) in current folder
vcf_files <- list.files(pattern = "\\.vcf\\.gz$")

for (vcf_file in vcf_files) {
  cat("Processing file:", vcf_file, "\n")

  # Read VCF
  vcf <- tryCatch(read.vcfR(vcf_file, verbose = FALSE), error = function(e) NULL)
  if (is.null(vcf)) next
  
  # Keep biallelic variants only
  vcf <- vcf[is.biallelic(vcf), ]
  if (nrow(vcf) == 0) next
  
  # Extract genotype data and recode
  gt <- extract.gt(vcf)
  gt <- as.data.frame(gt, stringsAsFactors = FALSE)
  gt[] <- lapply(gt, function(x) {
    x <- gsub("^0/0$", "1", x)
    x <- gsub("^1/1$", "2", x)
    x <- gsub("^(0/1|1/0)$", "3", x)
    as.numeric(replace(x, is.na(x) | x == ".", 0))
  })
  
  # Match sample IDs to phenotype info
  sample_ids <- colnames(gt)
  sample_info <- merge(data.frame(IDs = sample_ids), FullSpeciesInfo, by = "IDs", all.x = TRUE, sort = FALSE)
  sample_info$SexShort <- ifelse(sample_info$Sex == "Male", "M",
                                 ifelse(sample_info$Sex == "Female", "F", "U"))
  sample_info$name <- paste(sample_info$SexShort, sample_info$SampleName, sep = "-")
  colnames(gt) <- sample_info$name
  
  hetmatrix <- as.matrix(gt)
  
  # Identify male and female columns
  female_cols <- grep("^F-", colnames(hetmatrix))
  male_cols <- grep("^M-", colnames(hetmatrix))
  
  total_snps <- nrow(hetmatrix)
  count_fhet_mhom <- 0
  count_mhet_fhom <- 0
  
  # Count SNPs by sex-biased heterozygosity
  for (i in seq_len(total_snps)) {
    female_vals <- hetmatrix[i, female_cols]
    male_vals <- hetmatrix[i, male_cols]
    
    if (mean(female_vals == 3, na.rm = TRUE) >= 0.5 && mean(male_vals %in% c(1, 2), na.rm = TRUE) >= 0.5) {
      count_fhet_mhom <- count_fhet_mhom + 1
    }
    if (mean(male_vals == 3, na.rm = TRUE) >= 0.5 && mean(female_vals %in% c(1, 2), na.rm = TRUE) >= 0.5) {
      count_mhet_fhom <- count_mhet_fhom + 1
    }
  }
  
  # Extract chromosome, start, stop from filename (assumes format: Chr10_100_200.vcf.gz)
  vcf_info <- strsplit(basename(vcf_file), "[_.]")[[1]]
  chromosome <- gsub("Chr", "", vcf_info[1])
  start <- vcf_info[2]
  stop <- vcf_info[3]
  
  # Write summary counts
  summary_file <- paste0(Species, "_Chr", chromosome, "_", start, "_", stop, "_SexSample_het_hom_counts.txt")
  writeLines(
    c(
      paste("File:", vcf_file),
      paste("Total SNPs:", total_snps),
      paste("≥50% Female HET & ≥50% Male HOM:", count_fhet_mhom),
      paste("≥50% Male HET & ≥50% Female HOM:", count_mhet_fhom)
    ),
    con = summary_file
  )
  
  # Column annotation for heatmap
  col_annotation <- data.frame(Sex = sample_info$Sex)
  rownames(col_annotation) <- sample_info$name
  
  # Output files in all formats
  output_files <- list(
    pdf = paste0("pdf/", Species, "_heatmap_Chr", chromosome, "_", start, "_", stop, "_SexSample.pdf"),
    png = paste0("png/", Species, "_heatmap_Chr", chromosome, "_", start, "_", stop, "_SexSample.png"),
    tiff = paste0("tiff/", Species, "_heatmap_Chr", chromosome, "_", start, "_", stop, "_SexSample.tiff")
  )
  
  for (fmt in names(output_files)) {
    if (fmt == "pdf") pdf(output_files[[fmt]], width = 8.27, height = 11.69)
    if (fmt == "png") png(output_files[[fmt]], width = 8.27 * 300, height = 11.69 * 300, res = 300)
    if (fmt == "tiff") tiff(output_files[[fmt]], width = 8.27 * 300, height = 11.69 * 300, res = 300, compression = "lzw")
    
    pheatmap(
      hetmatrix,
      color = my_palette,
      cluster_rows = FALSE,
      show_rownames = TRUE,
      annotation_col = col_annotation,
      annotation_colors = annotation_colors,
      fontsize = 8,
      fontsize_col = 12,
      cellwidth = 14,
      legend_breaks = c(0, 1, 2, 3),
      legend_labels = c("NA", "homref", "homalt", "het"),
      main = paste("Genotype heatmap of", Species_name, "Chr", chromosome, ":", start, "-", stop)
    )
    dev.off()
  }
}
