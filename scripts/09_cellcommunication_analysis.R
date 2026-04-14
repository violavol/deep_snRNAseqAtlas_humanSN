library(Seurat)
library(CellChat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tibble)

source("R/cellchat_utils.R")

# Input files
atlas_rdata_path <- "data/processed/sn_atlas_annotated_subtype.RData"
dan_markers_path <- "Wtest_DaNs"
pd_risk_path <- "MAGMA_pd_risk_odcs2"

# Output directory
output_dir <- "results/cellcomm/"
#dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Analysis settings
disease_col <- "Disease"
celltype_col <- "CellSubType"
ctrl_label <- "CTR"
pd_label <- "PD_B5-6"

downsample_n <- 150
min_cells <- 10
n_pcs <- 10
cluster_resolution <- 0.6

# Enrichment test settings
n_top_dan_genes <- 1000
n_top_pd_risk_genes <- 200
bg_start <- 201
bg_end <- 1800
n_perm <- 10000

# Example pair
source_cell <- "ODC_2"
target_cell <- "DaN_3"

load(atlas_rdata_path)

cellchat_ctr <- prepare_cellchat_object(
  seurat_obj = sn_combined,
  disease_value = ctrl_label,
  disease_col = disease_col,
  celltype_col = celltype_col,
  downsample_n = downsample_n
)

cellchat_pd <- prepare_cellchat_object(
  seurat_obj = sn_combined,
  disease_value = pd_label,
  disease_col = disease_col,
  celltype_col = celltype_col,
  downsample_n = downsample_n
)

cellchat_ctr <- run_cellchat_pipeline(cellchat_ctr, min_cells = min_cells)
cellchat_pd  <- run_cellchat_pipeline(cellchat_pd,  min_cells = min_cells)

ctr_tables <- export_cellchat_tables(
  cellchat_obj = cellchat_ctr,
  prefix = "SNatlas_CTR_downsampled",
  output_dir = output_dir
)

pd_tables <- export_cellchat_tables(
  cellchat_obj = cellchat_pd,
  prefix = "SNatlas_PD_downsampled",
  output_dir = output_dir
)

df_net_ctr <- ctr_tables$df_net
df_net_pd  <- pd_tables$df_net

plot_signaling_role_heatmaps(cellchat_ctr)
plot_signaling_role_heatmaps(cellchat_pd)

cellchat_merged <- merge_cellchat_objects(cellchat_ctr, cellchat_pd)

p_count <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2))
p_weight <- compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "weight"
)

p_count + p_weight

p_heat_count <- netVisual_heatmap(cellchat_merged)
p_heat_weight <- netVisual_heatmap(cellchat_merged, measure = "weight")

p_heat_count + p_heat_weight

save(cellchat_merged,ctr_tables, file="results/cellcomm/cellchat_res.RData")

# Enrichment test for PD genetic risk-associated genes

wtest_dans <- read.delim(dan_markers_path, header = TRUE)
pd_risk_odc2 <- read.delim(pd_risk_path, header = TRUE)

dan3_ranked <- wtest_dans[order(-wtest_dans$DaN_3), , drop = FALSE]

top_dan_genes <- rownames(dan3_ranked)[seq_len(min(n_top_dan_genes, nrow(dan3_ranked)))]
top_pd_risk_genes <- pd_risk_odc2$ID[seq_len(min(n_top_pd_risk_genes, nrow(pd_risk_odc2)))]

background_pd_risk_genes <- pd_risk_odc2$ID[bg_start:min(bg_end, nrow(pd_risk_odc2))]

enrichment_result <- run_interaction_enrichment_test(
  cellcomm = df_net_ctr,
  dan_genes = top_dan_genes,
  pd_risk_genes = background_pd_risk_genes,
  source_cell = source_cell,
  target_cell = target_cell,
  n_perm = n_perm
)

enrichment_result
