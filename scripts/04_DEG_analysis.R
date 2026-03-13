library(Seurat)
library(dplyr)
library(here)

# Load subclustered object
load(here("data/processed/DA_subclusters.RData"))

# Differential expression
deg_results <- list()
for (clust in levels(Idents(DA_cells))) {
  deg_results[[clust]] <- FindMarkers(DA_cells, ident.1 = clust, min.pct = 0.25)
}

# Save DEG results
save(deg_results, file = here("results/DEG_tables/DA_DEG_results.RData"))
