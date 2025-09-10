# Load necessary libraries
library(ggtree)
library(treeio)
library(ggplot2)
library(stringr)
library(phangorn)

# Set base directory
base_dir <- file.path("Candidate_regions_iqtree", Species)  # Species variable defined below

# Create output directory
output_dir <- file.path(base_dir, "Collapsed_Tree_Plots")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Define species (change this to process different species)
Species <- "RhinePhoxinus"            # folder name
Species_name <- "Phoxinus phoxinus"  # species name for plot titles
Species_color <- "#66A61E"            # color for plots

# Load phenotype info (no header assumed)
pheno_file <- file.path(base_dir, paste0(Species, ".pheno.txt"))
FullSpeciesInfo <- read.delim(pheno_file, header = FALSE, stringsAsFactors = FALSE)
colnames(FullSpeciesInfo) <- c("IDs", "Species", "Sex", "Drainage", "SampleName")

# Convert Sex to M/F for consistency
FullSpeciesInfo$Sex <- ifelse(toupper(FullSpeciesInfo$Sex) == "FEMALE", "F",
                              ifelse(toupper(FullSpeciesInfo$Sex) == "MALE", "M", NA))

# Define color mapping for sexes
sex_color_mapping <- c(M = "#4D83C5", F = "#D62221", Unknown = "grey50")

# Function to plot tree
plot_tree <- function(tree_file, root_by_midpoint = TRUE) {
  file_name <- basename(tree_file)
  match <- str_match(file_name, "tree_(.+)_Chr(\\d+)_(\\d+)_(\\d+)\\.treefile")
  
  if (any(is.na(match))) {
    warning(paste("Skipping file with unmatched pattern:", file_name))
    return(NULL)
  }
  
  species <- match[2]
  chr <- match[3]
  start <- match[4]
  stop <- match[5]
  
  phylo_tree <- tryCatch(read.tree(tree_file), error = function(e) {
    message(paste("Error reading tree file:", tree_file, "\n", e))
    return(NULL)
  })
  if (is.null(phylo_tree)) return(NULL)
  
  if (root_by_midpoint) {
    phylo_tree <- phangorn::midpoint(phylo_tree)
  }
  
  polytree <- as.polytomy(phylo_tree, feature = "node.label", fun = function(x) as.numeric(x) < 50)
  
  new_labels <- sapply(polytree$tip.label, function(full_name) {
    parts <- strsplit(full_name, "_")[[1]]
    if (length(parts) < 5) return(full_name)
    
    sex <- ifelse(tolower(parts[1]) == "female", "F",
                  ifelse(tolower(parts[1]) == "male", "M", "U"))
    id_to_match <- paste(parts[2], parts[3], sep = "_")
    idx <- which(FullSpeciesInfo$IDs == id_to_match)
    sample_name <- ifelse(length(idx) == 1, FullSpeciesInfo$SampleName[idx], id_to_match)
    hap_tag <- parts[4]
    paste(sex, sample_name, hap_tag, sep = "-")
  })
  
  polytree$tip.label <- new_labels
  tip_groups <- substr(new_labels, 1, 1)
  tip_groups[!tip_groups %in% c("M", "F")] <- "Unknown"
  
  tip_data <- data.frame(label = polytree$tip.label, group = tip_groups)
  
  tree_plot <- ggtree(polytree, layout = "rectangular") %<+% tip_data +
    geom_tiplab(aes(label = label, color = group,
                    fontface = ifelse(group == "F", "italic", "plain")),
                size = 3, hjust = -0.3, vjust = 0.1) +
    scale_color_manual(values = sex_color_mapping) +
    geom_text2(aes(subset = !isTip, label = label), hjust = -0.3, size = 3) +
    geom_treescale(x = 0, y = -0.1, width = 0.04) +
    theme_tree() +
    labs(
      title = bquote("Haplotype tree for species: " ~ italic(.(Species_name))),
      subtitle = paste("Chr", chr, "Region:", start, "-", stop)
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(color = Species_color),
      plot.background = element_rect(fill = "white")
    )
  
  base_name <- paste0("tree_", species, "_Chr", chr, "_", start, "_", stop)
  ggsave(file.path(output_dir, paste0(base_name, ".png")), tree_plot, width = 14, height = 10, dpi = 300)
  ggsave(file.path(output_dir, paste0(base_name, ".pdf")), tree_plot, width = 14, height = 10)
  ggsave(file.path(output_dir, paste0(base_name, ".tiff")), tree_plot, width = 14, height = 10, dpi = 300)
}

# Run on all tree files in base_dir
tree_files <- list.files(base_dir, pattern = "\\.treefile$", full.names = TRUE, recursive = TRUE)
for (tree_file in tree_files) {
  tryCatch({
    plot_tree(tree_file, root_by_midpoint = TRUE)
    plot_tree(tree_file, root_by_midpoint = FALSE)
  }, error = function(e) {
    message(paste("Error processing file:", tree_file, "\n", e))
  })
}
