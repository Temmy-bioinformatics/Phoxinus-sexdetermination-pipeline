# Clear workspace
rm(list = ls())

# Load required libraries
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(patchwork)

# --------- SET BASE WORKING DIRECTORY ----------
# Change this to your base directory where species folders exist
base_dir <- "/path/to/Larger_candidate_regions_pca"

# --------- SPECIES SETTINGS ----------
species_data <- list(
  list(
    name = "Phoxinus phoxinus",
    directory = "RhinePhoxinus",  # Folder with MDS distance and metadata files
    color = "#66A61E",
    mdist_files = c("PCA_3_47140000_47160000.mdist",
                    "PCA_12_21270000_21280000.mdist")  # MDS distance matrix files
  ),
  list(
    name = "Phoxinus csikii",
    directory = "DanubeCsikii",
    color = "#D95E02",
    mdist_files = c("PCA_3_1490000_1499081.mdist")
  )
)

# Custom colors and shapes for sex groups in plots
color_pairs <- c("#D62221", "#4D83C5")  # Female = maroon, Male = turquoise
shape_palette <- c(15, 0)  # Female = filled square, Male = open square

# --------- FUNCTION TO GENERATE MDS PLOTS ----------
generate_mds_plot <- function(species_info, single_plot = FALSE) {
  # Compose species directory path
  species_dir <- file.path(base_dir, species_info$directory)
  setwd(species_dir)
  
  # Load metadata (sample info: IID, population, sex)
  # Expected file format: tab-delimited with no header, 3 columns: IID Population Sex
  metadata_file <- paste0(species_info$directory, "_metadata.id")
  data <- read.table(metadata_file, header = FALSE)
  colnames(data) <- c("IID", "Population", "Sex")
  famInd <- data.frame(IID = data$IID, Population = data$Population, Sex = data$Sex)
  
  # Select MDS distance files to process
  mdist_files <- species_info$mdist_files
  if (single_plot) {
    mdist_files <- mdist_files[1]  # Only first file if single_plot=TRUE
  }
  
  plot_list <- list()
  
  for (i in seq_along(mdist_files)) {
    mdist_file <- mdist_files[i]
    
    # Extract chromosome and genomic region from filename
    file_parts <- unlist(strsplit(mdist_file, split = "_|\\."))
    Chr <- paste0("Chr", file_parts[2])
    Region <- gsub("_", "-", paste0(file_parts[3], "_", file_parts[4]))
    
    # Load MDS distance matrix file (no header, square matrix)
    dist_populations <- read.table(mdist_file, header = FALSE)
    
    # Perform classical MDS (k=5 dimensions)
    mds_populations <- cmdscale(dist_populations, eig = TRUE, k = 5)
    
    # Combine eigenvectors with metadata
    eigenvec_populations <- cbind(famInd, mds_populations$points)
    colnames(eigenvec_populations)[4:8] <- paste0("Dim", 1:5)
    
    # Percent variance explained by each dimension
    eigen_percent <- round((mds_populations$eig / sum(mds_populations$eig)) * 100, 2)
    
    # Make Sex a factor for consistent plotting
    eigenvec_populations$Sex <- factor(eigenvec_populations$Sex, levels = c("Female", "Male"))
    
    # Create MDS scatter plot (Dim1 vs Dim2)
    mds <- ggplot(data = eigenvec_populations) +
      geom_point(aes(x = Dim1, y = Dim2, color = Sex, shape = Sex), size = 4, alpha = 0.7) +
      scale_shape_manual(values = shape_palette) +
      scale_color_manual(values = color_pairs) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
      geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
      labs(
        title = bquote("MDS plot for: " ~ italic(.(species_info$name))),
        subtitle = paste(Chr, ":", Region),
        x = paste0("Dimension 1 (", eigen_percent[1], "%)"),
        y = paste0("Dimension 2 (", eigen_percent[2], "%)")
      ) +
      theme_classic(base_size = 15) +
      theme(
        plot.title = element_text(size = 16, face = "bold", color = species_info$color),
        plot.subtitle = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none"
      )
    
    # Save individual plot if requested (example: Phoxinus csikii single plot)
    if (single_plot && species_info$name == "Phoxinus csikii") {
      ggsave("D.png", plot = mds, width = 10, height = 7)
    }
    
    plot_list[[i]] <- mds
  }
  
  # Return combined plot with shared legend, or single plot
  if (!single_plot) {
    combined_plot <- wrap_plots(plot_list, ncol = 2, guides = "collect") + 
      plot_annotation(
        title = paste("Combined MDS Plots for", species_info$name),
        theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
      )
    return(combined_plot)
  }
  
  return(plot_list[[1]])
}

# --------- GENERATE PLOTS ----------
phoxinus_phoxinus_combined <- generate_mds_plot(species_data[[1]])
phoxinus_csikii_plot <- generate_mds_plot(species_data[[2]], single_plot = TRUE)

# Extract legend from combined plot
legend <- get_legend(phoxinus_phoxinus_combined + theme(legend.position = "bottom"))

# Arrange plots without legend
combined_without_legend <- plot_grid(
  phoxinus_phoxinus_combined[[1]],
  phoxinus_phoxinus_combined[[2]],
  phoxinus_phoxinus_combined[[3]],
  phoxinus_csikii_plot,
  labels = c("A.", "B.", "C.", "D."),
  ncol = 2
)

# Add legend below plots
final_combined_plot <- plot_grid(
  combined_without_legend,
  legend,
  ncol = 1,
  rel_heights = c(0.9, 0.1)
)

# Save the final combined plot
ggsave("Phoxinus_MDS_Final_Combined.png", final_combined_plot, width = 16, height = 12, dpi = 300)

# Display the plot
print(final_combined_plot)
