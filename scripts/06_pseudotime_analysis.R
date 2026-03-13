library(Seurat)
library(monocle3)
library(here)

# Load subclustered dopaminergic neurons
load(here("data/processed/DA_subclusters.RData"))

# Convert Seurat object to Monocle3 CDS
cds <- as.cell_data_set(DA_cells)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

# Order cells along pseudotime
cds <- order_cells(cds)
DA_cells$pseudotime <- pseudotime(cds)

# Save updated Seurat object with pseudotime
save(DA_cells, file = here("data/processed/DA_subclusters_pseudotime.RData"))
