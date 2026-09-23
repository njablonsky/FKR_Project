library(ggplot2)
library(dplyr)

# 1. Load your data
df <- read.csv("marker_location.txt", sep="\t", header=TRUE)

# 2. Maintain Input Order
genome_order <- unique(df$Genome)
df$Genome <- factor(df$Genome, levels = rev(genome_order))

# 3. CRITICAL: Sort by Subfam, then by the visual rank of the Genome
# We convert Genome to numeric (which follows the rev(genome_order) levels)
# and sort descending so it follows Top -> Bottom
df <- df %>%
  mutate(genome_rank = as.numeric(Genome)) %>%
  arrange(Subfam, desc(genome_rank))

# 4. Prepare Backbone Data
backbone_df <- data.frame(
  Genome = rep(genome_order, each = 2),
  x = rep(c(0, 100), length(genome_order))
)
backbone_df$Genome <- factor(backbone_df$Genome, levels = rev(genome_order))

# 5. Create the Plot
ggplot() +
  # Backbones (Horizontal lines 0-100)
  geom_line(data = backbone_df, aes(x = x, y = Genome, group = Genome), 
            color = "grey90", linewidth = 0.5) +
  
  # USE geom_path INSTEAD OF geom_line
  # geom_path follows the row order of the dataframe strictly
  geom_path(data = df, aes(x = Location, y = Genome, group = Subfam, color = Subfam), 
            alpha = 0.4, linewidth = 0.4) +
  
  # Subfam points
  geom_point(data = df, aes(x = Location, y = Genome, color = Subfam), size = 1.5) +
  
  # Labels at the top
  # geom_text(data = subset(df, Genome == genome_order[1]), 
            # aes(x = Location, y = Genome, label = Subfam), 
            # vjust = -1.5, angle = 45, hjust = 0, 
            # size = 4, fontface = "bold") +
  
  # Axis and Scale
  scale_x_continuous(limits = c(0, 100), expand = c(0.01, 0.05)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 50, b = 10, l = 10)
  ) +
  coord_cartesian(clip = "off") + 
  labs(x = "Genomic location normalized by total protein counts", y = NULL)

## Plot gene order conservation

library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# 1. Load your pairwise data
data <- read.csv("pairwise_conservation_stats.csv")

# 2. Get unique list of all genomes to define the full matrix space
all_genomes <- sort(unique(c(data(Genome_A), data(Genome_B)))

# 3. Create the Self-Comparison (Diagonal) Data
self_comp <- data.frame(
  Genome_A = all_genomes,
  Genome_B = all_genomes,
  Shared_Markers = NA, 
  Kendall_Tau = 1.0,
  P_Value = 0.0
)

# 4. Combine: Original + Flipped + Self
data_sym <- bind_rows(
  data,                                              
  data %>% rename(Genome_A = Genome_B, Genome_B = Genome_A),
  self_comp                                          
) %>% 
  distinct(Genome_A, Genome_B, .keep_all = TRUE)

# 5. Pivot and FORCE factor levels to keep Row/Col order identical
tau_matrix <- data_sym %>%
  mutate(Genome_A = factor(Genome_A, levels = all_genomes),
         Genome_B = factor(Genome_B, levels = all_genomes)) %>%
  select(Genome_A, Genome_B, Kendall_Tau) %>%
  pivot_wider(names_from = Genome_B, values_from = Kendall_Tau, names_sort = TRUE) %>%
  tibble::column_to_rownames("Genome_A") %>%
  as.matrix()

p_matrix <- data_sym %>%
  mutate(Genome_A = factor(Genome_A, levels = all_genomes),
         Genome_B = factor(Genome_B, levels = all_genomes)) %>%
  select(Genome_A, Genome_B, P_Value) %>%
  pivot_wider(names_from = Genome_B, values_from = P_Value, names_sort = TRUE) %>%
  tibble::column_to_rownames("Genome_A") %>%
  as.matrix()

# 6. Final Polish
tau_matrix[is.na(tau_matrix)] <- 0
p_matrix[is.na(p_matrix)] <- 1

# 7. Define Color Scale (Standard: Red for Inversion, Blue for Synteny)
col_fun = colorRamp2(c(-1, 0, 1), c("#4575B4", "white", "#D73027"))

# 8. Plotting
pdf("Genomic_Synteny_Heatmap_Fixed.pdf", width = 16, height = 14)

ht=Heatmap(tau_matrix, 
        name = "Kendall Tau", 
        col = col_fun,
        column_title = "Gene Order Conservation (Significance marked with *)",
        row_names_gp = gpar(fontsize = 6), # Small font for large genome lists
        column_names_gp = gpar(fontsize = 6),
        show_column_names = TRUE,
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "complete",
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "complete",
        # Logic for concise labeling
        cell_fun = function(j, i, x, y, width, height, fill) {
          p_val <- p_matrix[i, j]
          # Only label with a single star if p < 0.05
          if (!is.na(p_val) && p_val < 0.05) {
            grid.text("*", x, y, gp = gpar(fontsize = 8, col = "black"))
          }
        },
        heatmap_legend_param = list(
          title = "Tau Score",
          at = c(-1, 0, 1),
          labels = c("Inverted", "Random", "Conserved")
        ))

# Use draw() to explicitly control the legend position and padding
draw(ht, 
     heatmap_legend_side = "right", 
     # padding: c(bottom, left, top, right) 
     # Adding 20 units to the right moves the legend away from the labels
     padding = unit(c(0, 20, 20, 50), "mm"))

dev.off()

ht

# Save as a comma-separated values file
write.csv(data_sym, "45_subfam_location_conservation.csv", row.names = FALSE)
