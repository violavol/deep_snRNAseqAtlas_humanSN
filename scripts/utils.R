library(ggplot2)
library(dplyr)
library(Seurat)

# Calculate Jaccard index between two gene sets
jaccard <- function(a, b) {
  inter <- length(intersect(a, b))
  union <- length(a) + length(b) - inter
  inter / union
}

# Save plot utility
save_plot <- function(plot_obj, filename, width = 6, height = 4) {
  ggsave(filename, plot = plot_obj, width = width, height = height)
}

# DotPlot wrapper
dotplot_rot <- function(seurat_obj, features, cluster_ids = NULL, scale = 8) {
  DotPlot(seurat_obj, features = features, idents = cluster_ids, dot.scale = scale) +
    RotatedAxis() + theme_minimal()
}
