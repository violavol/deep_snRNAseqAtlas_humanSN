library(Seurat)
library(dplyr)
library(here)

# Load raw data
raw_file <- here("data/raw/sn_atlas_unprocessed.RData")
if (!file.exists(raw_file)) stop("Raw data file missing!")

load(raw_file) 
sn_atlas <- sn_atlas_unprocessed
rm(sn_atlas_unprocessed)

# Split by disease and normalize
sn_list <- SplitObject(sn_atlas, split.by = "Disease")
sn_list <- lapply(sn_list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Integrate
features <- SelectIntegrationFeatures(sn_list)
sn_list <- lapply(sn_list, function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
})

anchors <- FindIntegrationAnchors(sn_list, anchor.features = features, reduction = "rpca", k.anchor = 20)
sn_combined <- IntegrateData(anchorset = anchors)

DefaultAssay(sn_combined) <- "integrated"
sn_combined <- ScaleData(sn_combined)
sn_combined <- RunPCA(sn_combined, npcs = 40)
sn_combined <- RunUMAP(sn_combined, reduction = "pca", dims = 1:40)
sn_combined <- FindNeighbors(sn_combined, reduction = "pca", dims = 1:40)
sn_combined <- FindClusters(sn_combined, resolution = 0.5)

# Save processed object
save(sn_combined, file = here("data/processed/sn_atlas_integrated.RData"))
