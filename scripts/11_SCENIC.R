# 11_SCENIC.R
# SCENIC transcription-factor regulon analysis on DaN (CTR only).
# Project: deep_snRNAseqAtlas_humanSN — Volpato et al. (2026)

suppressPackageStartupMessages({
  library(SCENIC)
  library(AUCell)
  library(RcisTarget)
  library(SCopeLoomR)
  library(KernSmooth)
  library(BiocParallel)
  library(ggplot2)
  library(data.table)
  library(grid)
  library(ComplexHeatmap)
  library(Seurat)
  library(doRNG)
  library(here)
})
source(here("scripts/utils.R"))

set.seed(42)

# ── Parameters ────────────────────────────────────────────────────────────────
CTRL_LABEL  <- "CTR"
N_CORES     <- 10
SEED        <- 42
DB_DIR      <- here("data/external/cisTarget_databases")  # portable path

# ── Load data ─────────────────────────────────────────────────────────────────
load(here("data/processed/sn_atlas_annotated_subtype.RData"))

sn_dan_ctr <- subset(sn_combined,
                     subset = CellType == "DaN" & Disease == CTRL_LABEL)
message("DaN CTR cells for SCENIC: ", ncol(sn_dan_ctr))

genes_dan_ctr     <- as.matrix(GetAssayData(sn_dan_ctr, slot = "counts"))
metadata_dan_ctr  <- as.data.frame(sn_dan_ctr@meta.data)   # was wrongly "tmp" in original

# ── Initialise SCENIC ─────────────────────────────────────────────────────────
org <- "hgnc"
data(defaultDbNames)
dbs <- defaultDbNames[[org]]

scenicOptions <- initializeScenic(
  org          = org,
  dbDir        = DB_DIR,
  dbs          = dbs,
  datasetTitle = "DaN_CTR",
  nCores       = N_CORES
)
scenicOptions@inputDatasetInfo$cellInfo  <- metadata_dan_ctr
scenicOptions@inputDatasetInfo$colVars   <- metadata_dan_ctr$CellSubType

data(list = "motifAnnotations_hgnc_v9", package = "RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9

# ── Co-expression network ─────────────────────────────────────────────────────
set.seed(SEED)
runCorrelation(genes_dan_ctr, scenicOptions)

genes_dan_ctr_log <- log2(genes_dan_ctr + 1)
set.seed(SEED)
runGenie3(genes_dan_ctr_log, scenicOptions)
rm(genes_dan_ctr); gc()

# ── Build and score GRN ───────────────────────────────────────────────────────
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
saveRDS(scenicOptions, file = here("results/SCENIC/scenicOptions_1.Rds"))

scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
saveRDS(scenicOptions, file = here("results/SCENIC/scenicOptions_2.Rds"))

scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, genes_dan_ctr_log)
saveRDS(scenicOptions, file = here("results/SCENIC/scenicOptions_3.Rds"))

# ── Regulon activity per cell sub-type ───────────────────────────────────────
regulonAUC_dans <- loadInt(scenicOptions, "aucell_regulonAUC")

# Use sn_dan_ctr — NOT an undefined "tmp" object
Idents(sn_dan_ctr) <- sn_dan_ctr$CellSubType
cell_info <- data.frame(seuratCluster = Idents(sn_dan_ctr),
                        row.names     = colnames(sn_dan_ctr))

regulon_by_celltype <- sapply(
  split(rownames(cell_info), cell_info[, 1]),
  function(cells) rowMeans(getAUC(regulonAUC_dans)[, cells, drop = FALSE])
)

regulon_scaled <- t(scale(t(regulon_by_celltype), center = TRUE, scale = TRUE))

p_heat <- ComplexHeatmap::Heatmap(
  regulon_scaled,
  name          = "Regulon activity",
  row_names_gp  = grid::gpar(fontsize = 5)
)

pdf(here("figures/Figure_6/SCENIC_regulon_heatmap.pdf"), width = 8, height = 12)
draw(p_heat)
dev.off()

# ── Top regulators per cell type ──────────────────────────────────────────────
top_regulators <- reshape2::melt(regulon_scaled)
colnames(top_regulators) <- c("Regulon", "CellType", "RelativeActivity")
top_regulators <- top_regulators[top_regulators$RelativeActivity > 0, ]

# ── Save ──────────────────────────────────────────────────────────────────────
ensure_dirs(here("results/SCENIC"))
save(regulonAUC_dans, regulon_by_celltype, regulon_scaled, top_regulators,
     file = here("results/SCENIC/SCENIC_results.RData"))
message("SCENIC results saved.")

log_session("11_SCENIC")
