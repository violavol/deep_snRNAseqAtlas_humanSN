# 05_DTU_analysis.R
# Differential transcript usage (DTU) using fishpond/swish, plus
# isoform-level pseudotime analysis along the OPC → ODC_2 trajectory.
# Project: deep_snRNAseqAtlas_humanSN — Volpato et al. (2026)

library(tximeta)
library(fishpond)
library(scran)
library(biomaRt)
library(slingshot)
library(Seurat)
library(switchde)
library(here)
source(here("scripts/utils.R"))

set.seed(42)

# ── Parameters ────────────────────────────────────────────────────────────────
CTRL_LABEL  <- "CTR"
PD_LABEL    <- "PD_B5-6"
MIN_COUNT   <- 3
MIN_N       <- 10
N_PERMS     <- 64
N_DIMS_ISO  <- 10
CLUSTER_RES <- 0.6
MEAN_EXPR_THRESHOLD   <- 0.1
NONZERO_FRAC_THRESHOLD <- 0.2
SEED <- 42

# ── Load salmon quantifications ───────────────────────────────────────────────
coldata <- read.delim(here("data/processed/coldata_ALL_CTRandPD56_modified_nfs"),
                      header = TRUE, stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(SummarizedExperiment))
y <- tximeta(coldata, dropInfReps = TRUE)

# ── Differential transcript expression ───────────────────────────────────────
y <- labelKeep(y, minCount = MIN_COUNT, minN = MIN_N)
y <- y[mcols(y)$keep, ]

assays(y) <- lapply(assays(y), as.matrix)   # dense needed for scran
y <- scaleInfReps(y, lengthCorrect = FALSE, sfFun = computeSumFactors)
y <- swish(y, x = "condition", quiet = TRUE)

# ── Differential isoform usage ────────────────────────────────────────────────
iso     <- isoformProportions(y)
iso     <- swish(iso, x = "condition", nperms = N_PERMS)
dtu_df  <- as.data.frame(mcols(iso)[, c("log2FC", "qvalue", "gene", "tx_id")])

# ── OPC → ODC_2 pseudotime on isoform data ───────────────────────────────────
load(here("data/processed/pseudotime_isoforms_OPC_ODCs_CTRandPD56.RData"))

isoform_mat <- as.matrix(assays(y)[[2]])
meta_iso <- data.frame(
  celltype = y$cellSubType,
  disease  = y$condition,
  cell     = y$Barcode,
  sample   = y$names,
  row.names = paste(y$Barcode, y$names, sep = "_")
)
colnames(isoform_mat) <- rownames(meta_iso)

run_opc_odc_pseudotime <- function(isoform_mat, meta_iso, disease_val, seed = SEED) {
  meta_sub <- meta_iso[meta_iso$disease == disease_val, ]
  mat_sub  <- isoform_mat[, rownames(meta_sub), drop = FALSE]

  obj <- CreateSeuratObject(counts = mat_sub, min.cells = 0,
                            min.features = 0, meta.data = meta_sub)

  set.seed(seed)
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                              nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:N_DIMS_ISO, verbose = FALSE)

  set.seed(seed)
  obj <- FindClusters(obj, resolution = CLUSTER_RES, verbose = FALSE)

  set.seed(seed)
  obj <- RunUMAP(obj, dims = 1:N_DIMS_ISO, seed.use = seed, verbose = FALSE)

  sce <- as.SingleCellExperiment(obj, assay = "RNA")
  sce <- slingshot(sce, reducedDim = "UMAP", clusterLabels = "seurat_clusters")

  list(obj = obj, sce = sce)
}

res_ctr <- run_opc_odc_pseudotime(isoform_mat, meta_iso, CTRL_LABEL)
res_pd  <- run_opc_odc_pseudotime(isoform_mat, meta_iso, PD_LABEL)

meta_iso_ctr <- meta_iso[meta_iso$disease == CTRL_LABEL, ]
meta_iso_ctr$pseudotime <- res_ctr$sce$slingPseudotime_1

meta_iso_pd <- meta_iso[meta_iso$disease == PD_LABEL, ]
meta_iso_pd$pseudotime <- res_pd$sce$slingPseudotime_1

# ── SwitchDE along OPC→ODC_2 trajectory (CTR) ────────────────────────────────
iso_ctr_log <- log(as.matrix(GetAssayData(res_ctr$obj, slot = "counts")) + 1)
iso_ctr_log <- iso_ctr_log[
  rowMeans(iso_ctr_log) > MEAN_EXPR_THRESHOLD &
    rowMeans(iso_ctr_log > 0) > NONZERO_FRAC_THRESHOLD, ]

sde_iso_opc_odc2_ctr <- switchde(iso_ctr_log, res_ctr$sce$slingPseudotime_1)

# ── Save ──────────────────────────────────────────────────────────────────────
ensure_dirs(here("results/pseudotime"))

save(dtu_df, meta_iso_ctr, meta_iso_pd, sde_iso_opc_odc2_ctr,
     file = here("results/pseudotime/OPC_ODCs_DTU_ps.RData"))
message("DTU and isoform pseudotime results saved.")

log_session("05_DTU_analysis")
