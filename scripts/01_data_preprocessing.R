# 01_data_preprocessing.R
# Integration of raw snRNA-seq data across disease conditions.
# Project: deep_snRNAseqAtlas_humanSN — Volpato et al. (2026)

library(Seurat)
library(dplyr)
library(here)
source(here("scripts/utils.R"))

# ── Parameters ────────────────────────────────────────────────────────────────
N_FEATURES    <- 2000   # variable features per sample
N_PCS         <- 40     # PCs for UMAP / neighbours
K_ANCHOR      <- 20     # anchors for rpca integration
CLUSTER_RES   <- 0.5    # Leiden resolution
SEED          <- 42

# ── Load raw data ─────────────────────────────────────────────────────────────
raw_file <- here("data/raw/sn_atlas_unprocessed.RData")
if (!file.exists(raw_file)) stop("Raw data file missing: ", raw_file)

load(raw_file)
sn_atlas <- sn_atlas_unprocessed
rm(sn_atlas_unprocessed)
gc()

message("Loaded atlas: ", ncol(sn_atlas), " cells × ", nrow(sn_atlas), " genes")

# ── Per-sample normalisation and variable-feature selection ───────────────────
sn_list <- SplitObject(sn_atlas, split.by = "Disease")
rm(sn_atlas); gc()

sn_list <- lapply(sn_list, function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst",
                            nfeatures = N_FEATURES, verbose = FALSE)
  x
})

# ── rPCA integration ──────────────────────────────────────────────────────────
features <- SelectIntegrationFeatures(sn_list, nfeatures = N_FEATURES)

sn_list <- lapply(sn_list, function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
  x
})

anchors     <- FindIntegrationAnchors(sn_list,
                                      anchor.features = features,
                                      reduction       = "rpca",
                                      k.anchor        = K_ANCHOR)
sn_combined <- IntegrateData(anchorset = anchors)
rm(sn_list, anchors); gc()

# ── Dimensionality reduction & clustering ─────────────────────────────────────
DefaultAssay(sn_combined) <- "integrated"

set.seed(SEED)
sn_combined <- ScaleData(sn_combined, verbose = FALSE)
sn_combined <- RunPCA(sn_combined, npcs = N_PCS, verbose = FALSE)
sn_combined <- RunUMAP(sn_combined, reduction = "pca", dims = 1:N_PCS, seed.use = SEED)
sn_combined <- FindNeighbors(sn_combined, reduction = "pca", dims = 1:N_PCS)

set.seed(SEED)
sn_combined <- FindClusters(sn_combined, resolution = CLUSTER_RES)

message("Clustering complete: ", nlevels(Idents(sn_combined)), " clusters")

# ── Save ──────────────────────────────────────────────────────────────────────
ensure_dirs(here("data/processed"))
save(sn_combined, file = here("data/processed/sn_atlas_integrated.RData"))
message("Saved integrated object.")

log_session("01_data_preprocessing")
