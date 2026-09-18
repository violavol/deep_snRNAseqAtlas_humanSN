# 03_subtype_analysis.R
# Sub-clustering of dopaminergic neurons (DaN) and marker visualisation.
# The same integration recipe can be applied to other cell types by adjusting
# the CELL_TYPE parameter and, if needed, K_WEIGHT.
# Project: deep_snRNAseqAtlas_humanSN — Volpato et al. (2026)

library(Seurat)
library(here)
source(here("scripts/utils.R"))

# ── Parameters ────────────────────────────────────────────────────────────────
CELL_TYPE   <- "DaN"
N_FEATURES  <- 2000
N_PCS       <- 40
CLUSTER_RES <- 0.2      # lower resolution → broader sub-clusters
K_ANCHOR    <- 40
SEED        <- 42
OUTLIER_SAMPLE <- "14_133"   # excluded: outlier in UMAP and pseudotime

# DaN sub-type marker panel
DAN_MARKERS <- c(
  "SOX6", "GRIA3", "DCX", "AGTR1", "LMX1B", "GRIK3",
  "RET", "GFRA2", "PITX3", "CHRNA4", "SLC18A2", "SLC6A3",
  "TH", "KCNJ6", "ALDH1A1", "TMEM255A", "LGI1", "TMEFF2"
)

# ── Load annotated object ─────────────────────────────────────────────────────
load(here("data/processed/sn_atlas_annotated.RData"))

# ── Subset target cell type, remove outlier sample ───────────────────────────
da_cells <- subset(sn_combined, idents = CELL_TYPE)
da_cells <- subset(da_cells, subset = Sample_v2 != OUTLIER_SAMPLE)

message(sprintf("DaN cells after QC filter: %d", ncol(da_cells)))

# ── rPCA integration within DaN ───────────────────────────────────────────────
dan_list <- SplitObject(da_cells, split.by = "Disease")

dan_list <- lapply(dan_list, function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst",
                            nfeatures = N_FEATURES, verbose = FALSE)
  x
})

features <- SelectIntegrationFeatures(dan_list, nfeatures = N_FEATURES)

dan_list <- lapply(dan_list, function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
  x
})

anchors      <- FindIntegrationAnchors(dan_list,
                                       anchor.features = features,
                                       reduction       = "rpca",
                                       k.anchor        = K_ANCHOR)
dan_combined <- IntegrateData(anchorset = anchors)
rm(da_cells, dan_list, anchors); gc()

# ── Dimensionality reduction & sub-clustering ─────────────────────────────────
DefaultAssay(dan_combined) <- "integrated"

set.seed(SEED)
dan_combined <- ScaleData(dan_combined, verbose = FALSE)
dan_combined <- RunPCA(dan_combined, npcs = N_PCS, verbose = FALSE)
dan_combined <- RunUMAP(dan_combined, reduction = "pca", dims = 1:N_PCS, seed.use = SEED)
dan_combined <- FindNeighbors(dan_combined, reduction = "pca", dims = 1:N_PCS)

set.seed(SEED)
dan_combined <- FindClusters(dan_combined, resolution = CLUSTER_RES)

message("DaN sub-clusters: ", nlevels(Idents(dan_combined)))

# ── Marker dot-plot (clusters 0–3; 4–5 too small / disease-restricted) ───────
DefaultAssay(dan_combined) <- "RNA"

p_dot <- dotplot_rot(dan_combined, features = DAN_MARKERS,
                     cluster_ids = as.character(0:3))
save_plot(p_dot, here("figures/Figure_2/dotplot_DaN_markers.pdf"), width = 12, height = 5)

p_umap <- DimPlot(dan_combined, reduction = "umap",
                  split.by = "Disease", group.by = "seurat_clusters")
save_plot(p_umap, here("figures/Figure_2/umap_DaN_subclusters.pdf"), width = 16, height = 5)

# ── Transfer sub-cluster labels to the full atlas object ─────────────────────
# Use %in% for set membership — NOT == — to handle vector-length mismatches
sub_meta <- as.data.frame(dan_combined@meta.data)

# Initialise column (NA for non-DaN cells)
sn_combined$CellSubType <- NA_character_

for (cl in as.character(sort(unique(sub_meta$seurat_clusters)))) {
  cells_in_cluster <- rownames(sub_meta)[sub_meta$seurat_clusters == cl]
  sn_combined$CellSubType[colnames(sn_combined) %in% cells_in_cluster] <-
    paste0("DaN_", cl)
}

message("CellSubType distribution (DaN cells):")
print(table(sn_combined$CellSubType, useNA = "ifany"))

# ── Save ──────────────────────────────────────────────────────────────────────
save(sn_combined, dan_combined,
     file = here("data/processed/sn_atlas_annotated_subtype.RData"))
message("Saved subtype-annotated object.")

log_session("03_subtype_analysis")
