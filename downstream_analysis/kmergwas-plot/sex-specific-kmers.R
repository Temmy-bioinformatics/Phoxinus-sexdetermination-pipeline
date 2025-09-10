# Clear workspace
rm(list = ls())

# Load libraries
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpattern)

# Set working directory relative to project root
setwd("kmerGWAS_male_female_specific")

# Species metadata: names, directories, and colors
species_data <- list(
  list(name = "Phoxinus phoxinus", directory = "RhinePhoxinus", color = "#66A61E"),
  list(name = "Phoxinus csikii", directory = "DanubeCsikii", color = "#D95E02")
)

# Define sex‑specific aesthetics
sex_colors <- c("Male-specific" = "#4D83C5", "Female-specific" = "#D62221")
pattern_types <- c("Male-specific" = "none", "Female-specific" = "stripe")

# Initialize combined data frame
combined_counts <- data.frame()

# Load and process each species' k‑mer data
for (species_info in species_data) {
  species <- species_info$directory
  species_name <- species_info$name
  
  file_path <- file.path(species, "assoc_kmers_presab.txt")
  data <- read.table(file_path, header = TRUE, sep = "\t")
  
  # Identify columns by sex
  male_columns <- grep("Male", colnames(data), value = TRUE)
  female_columns <- grep("Female", colnames(data), value = TRUE)
  
  # Sex‑specific k‑mers: present in all of one sex, absent in all of the other
  male_kmers <- data[rowSums(data[male_columns] == 1) > 0 & rowSums(data[female_columns] == 0) == length(female_columns), ]
  female_kmers <- data[rowSums(data[female_columns] == 1) > 0 & rowSums(data[male_columns] == 0) == length(male_columns), ]
  
  # Append counts
  combined_counts <- rbind(combined_counts, data.frame(
    Species = species_name,
    Count = c(nrow(male_kmers), nrow(female_kmers)),
    Sex = c("Male-specific", "Female-specific")
  ))
}

# Factor ordering for plotting
combined_counts$Species <- factor(combined_counts$Species, levels = c("Phoxinus phoxinus", "Phoxinus csikii"))
combined_counts$Sex <- factor(combined_counts$Sex, levels = c("Male-specific", "Female-specific"))

# Function to create individual panel plots
generate_panel_plot <- function(data, species_name, species_color, panel_label, show_y_axis = TRUE) {
  ggplot(data, aes(x = Sex, y = Count, fill = Sex, pattern = Sex)) +
    geom_bar_pattern(
      stat = "identity",
      width = 0.3,
      pattern_fill = "black",
      pattern_colour = "black",
      pattern_density = 0.001,
      pattern_spacing = 0.05,
      pattern_angle = 45
    ) +
    scale_fill_manual(values = c("Male-specific" = species_color, "Female-specific" = species_color)) +
    scale_pattern_manual(values = pattern_types) +
    labs(
      title = bquote(.(panel_label) ~ "Sex-specific k-mers in" ~ italic(.(species_name))),
      x = NULL,
      y = if (show_y_axis) "Count of Sex-Specific k-mers" else NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0, size = 16, color = species_color),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y = if (show_y_axis) element_text(size = 13) else element_blank(),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# Generate plots per species
plot_phoxinus <- generate_panel_plot(
  combined_counts %>% filter(Species == "Phoxinus phoxinus"),
  "Phoxinus phoxinus", "#66A61E", "A.", TRUE
)
plot_csikii <- generate_panel_plot(
  combined_counts %>% filter(Species == "Phoxinus csikii"),
  "Phoxinus csikii", "#D95E02", "B.", FALSE
)

# Save individual plots
ggsave("Phoxinus_phoxinus_sex_specific_kmers.pdf", plot_phoxinus, width = 6, height = 4, dpi = 300)
ggsave("Phoxinus_phoxinus_sex_specific_kmers.png", plot_phoxinus, width = 10, height = 5, dpi = 600)
ggsave("Phoxinus_csikii_sex_specific_kmers.pdf", plot_csikii, width = 6, height = 5, dpi = 300)
ggsave("Phoxinus_csikii_sex_specific_kmers.png", plot_csikii, width = 10, height = 5, dpi = 600)

# Combine both species plots into one panel
final_plot <- plot_phoxinus + plot_csikii + plot_layout(ncol = 2)

# Save combined plot
ggsave("Phoxinus_sex_specific_kmers_bw_thinbars.pdf", final_plot, width = 12, height = 5, dpi = 300)
ggsave("Phoxinus_sex_specific_kmers_bw_thinbars.png", final_plot, width = 12, height = 5, dpi = 600)

# Display in R session
print(final_plot)
