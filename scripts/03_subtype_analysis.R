library(Seurat)
library(here)

# Load annotated object
load(here("data/processed/sn_atlas_annotated.RData"))

# Sub-cluster dopaminergic neurons
DA_cells <- subset(sn_combined, idents = "Dopaminergic")
DA_cells <- FindVariableFeatures(DA_cells)
DA_cells <- ScaleData(DA_cells)
DA_cells <- RunPCA(DA_cells)
DA_cells <- RunUMAP(DA_cells, dims = 1:20)
DA_cells <- FindNeighbors(DA_cells, dims = 1:20)
DA_cells <- FindClusters(DA_cells, resolution = 0.3)

# Save subclustered object
save(DA_cells, file = here("data/processed/DA_subclusters.RData"))
