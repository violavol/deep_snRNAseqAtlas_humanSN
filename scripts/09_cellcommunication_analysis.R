# 09_cellcommunication_analysis.R
# CellChat-based cell-cell communication analysis (CTR vs PD_B5-6) plus
# permutation enrichment test for PD genetic-risk genes on ODC_2 → DaN_3 interactions.
# Project: deep_snRNAseqAtlas_humanSN — Volpato et al. (2026)

library(Seurat)
library(CellChat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(here)
source(here("scripts/utils.R"))
source(here("scripts/cellchat_utils.R"))

# ── Parameters ────────────────────────────────────────────────────────────────
CTRL_LABEL        <- "CTR"
PD_LABEL          <- "PD_B5-6"
DISEASE_COL       <- "Disease"
CELLTYPE_COL      <- "CellSubType"
DOWNSAMPLE_N      <- 150
MIN_CELLS         <- 10
N_TOP_DAN_GENES   <- 1000
N_TOP_PD_RISK     <- 200
BG_START          <- 201
BG_END            <- 1800
N_PERM            <- 10000
SOURCE_CELL       <- "ODC_2"
TARGET_CELL       <- "DaN_3"

# ── Input / output paths ──────────────────────────────────────────────────────
OUTPUT_DIR  <- here("results/cellcomm")
ensure_dirs(OUTPUT_DIR)

# ── Load atlas ────────────────────────────────────────────────────────────────
load(here("data/processed/sn_atlas_annotated_subtype.RData"))

# ── Build CellChat objects ────────────────────────────────────────────────────
cellchat_ctr <- prepare_cellchat_object(
  seurat_obj   = sn_combined,
  disease_value = CTRL_LABEL,
  disease_col  = DISEASE_COL,
  celltype_col = CELLTYPE_COL,
  downsample_n = DOWNSAMPLE_N
)
cellchat_pd <- prepare_cellchat_object(
  seurat_obj   = sn_combined,
  disease_value = PD_LABEL,
  disease_col  = DISEASE_COL,
  celltype_col = CELLTYPE_COL,
  downsample_n = DOWNSAMPLE_N
)

# ── Run pipelines ─────────────────────────────────────────────────────────────
cellchat_ctr <- run_cellchat_pipeline(cellchat_ctr, min_cells = MIN_CELLS)
cellchat_pd  <- run_cellchat_pipeline(cellchat_pd,  min_cells = MIN_CELLS)

# ── Export tables ─────────────────────────────────────────────────────────────
ctr_tables <- export_cellchat_tables(cellchat_ctr,
                                     prefix     = "SNatlas_CTR_downsampled",
                                     output_dir = OUTPUT_DIR)
pd_tables  <- export_cellchat_tables(cellchat_pd,
                                     prefix     = "SNatlas_PD_downsampled",
                                     output_dir = OUTPUT_DIR)
df_net_ctr <- ctr_tables$df_net
df_net_pd  <- pd_tables$df_net

# ── Signalling role heatmaps ──────────────────────────────────────────────────
plot_signaling_role_heatmaps(cellchat_ctr)
plot_signaling_role_heatmaps(cellchat_pd)

# ── Merge & compare ───────────────────────────────────────────────────────────
cellchat_merged <- merge_cellchat_objects(cellchat_ctr, cellchat_pd)

p_count  <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2))
p_weight <- compareInteractions(cellchat_merged, show.legend = FALSE,
                                group = c(1, 2), measure = "weight")
p_count + p_weight

p_heat_count  <- netVisual_heatmap(cellchat_merged)
p_heat_weight <- netVisual_heatmap(cellchat_merged, measure = "weight")
p_heat_count + p_heat_weight

# ── PD risk enrichment test ───────────────────────────────────────────────────
wtest_dans <- read.delim(here("results/DEG_tables/Wtest_DaN_subtypes.tsv"),
                         header = TRUE)
pd_risk_odc2 <- read.delim(here("data/external/MAGMA_pd_risk_odcs2"),
                            header = TRUE)

dan3_ranked   <- wtest_dans[order(-wtest_dans$DaN_3), , drop = FALSE]
top_dan_genes <- rownames(dan3_ranked)[seq_len(min(N_TOP_DAN_GENES, nrow(dan3_ranked)))]

top_pd_risk_genes <- pd_risk_odc2$ID[seq_len(min(N_TOP_PD_RISK, nrow(pd_risk_odc2)))]
bg_pd_risk_genes  <- pd_risk_odc2$ID[BG_START:min(BG_END, nrow(pd_risk_odc2))]

enrichment_result <- run_interaction_enrichment_test(
  cellcomm    = df_net_ctr,
  dan_genes   = top_dan_genes,
  pd_risk_genes = bg_pd_risk_genes,
  source_cell = SOURCE_CELL,
  target_cell = TARGET_CELL,
  n_perm      = N_PERM
)
message(sprintf("Enrichment p-value: %.4f  (n_observed = %d)",
                enrichment_result$p_value, enrichment_result$n_observed))

# ── Save ──────────────────────────────────────────────────────────────────────
save(cellchat_merged, ctr_tables, pd_tables, enrichment_result,
     file = here("results/cellcomm/cellchat_res.RData"))
message("CellChat results saved.")

log_session("09_cellcommunication_analysis")
