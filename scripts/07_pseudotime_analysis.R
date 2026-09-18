# 07_pseudotime_analysis.R
# Pseudotime analysis for DaN (monocle3 + slingshot) and ODC (slingshot),
# plus DE genes along trajectory using switchDE and tradeSeq.
# Project: deep_snRNAseqAtlas_humanSN — Volpato et al. (2026)

library(Seurat)
library(monocle3)
library(slingshot)
library(tradeSeq)
library(switchde)
library(ggplot2)
library(reshape2)
library(igraph)
library(here)
source(here("scripts/utils.R"))

set.seed(42)

# ── Parameters ────────────────────────────────────────────────────────────────
SEED          <- 42
CTRL_LABEL    <- "CTR"
PD_LABEL      <- "PD_B5-6"
MONOCLE_DIMS  <- 3
MONOCLE_RES   <- 1e-5
PS_BREAKS     <- c(0, 5, 10, 15, 20)   # pseudotime interval cut-points
MEAN_EXPR_THR <- 0.1
NONZERO_THR   <- 0.2

# Helper: compute TPM-normalised, log-transformed expression
tpm_log <- function(counts, gene_lengths) {
  prot_cod <- gene_lengths[match(rownames(counts), gene_lengths[, 1]), ]
  mat <- sweep(counts, 1, STATS = prot_cod$Gene_Length / 1000, FUN = "/")
  mat[is.na(mat)] <- 0
  mat <- sweep(mat, 2, STATS = colSums(mat) / 1e6, FUN = "/")
  log(mat + 1)
}

# Helper: filter low-expressed genes before switchDE
filter_expr <- function(mat, mean_thr = MEAN_EXPR_THR, nonzero_thr = NONZERO_THR) {
  mat[rowMeans(mat) > mean_thr & rowMeans(mat > 0) > nonzero_thr, ]
}

# ── Load data ─────────────────────────────────────────────────────────────────
load(here("data/processed/sn_atlas_annotated_subtype.RData"))
gene_lengths <- read.delim(here("data/processed/gene_length"), header = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: DaN pseudotime with monocle3
# ══════════════════════════════════════════════════════════════════════════════

sn_dan <- subset(sn_combined, subset = CellType == "DaN")
message("DaN cells: ", ncol(sn_dan))

counts_dan  <- as.matrix(GetAssayData(sn_dan, slot = "counts"))
meta_dan    <- sn_dan@meta.data
gene_meta   <- data.frame(
  gene_short_name = rownames(counts_dan),
  genes           = rownames(counts_dan),
  row.names       = rownames(counts_dan)
)

cds_dan <- new_cell_data_set(counts_dan,
                             cell_metadata = meta_dan,
                             gene_metadata = gene_meta)
rm(counts_dan); gc()

set.seed(SEED)
cds_dan <- preprocess_cds(cds_dan, num_dim = MONOCLE_DIMS)
cds_dan <- reduce_dimension(cds_dan)
cds_dan <- cluster_cells(cds_dan, resolution = MONOCLE_RES)
cds_dan <- learn_graph(cds_dan)

plot_cells(cds_dan, color_cells_by = "Disease",
           label_branch_points = FALSE, label_roots = FALSE,
           label_leaves = FALSE, cell_size = 1)

set.seed(SEED)
cds_dan <- order_cells(cds_dan)

plot_cells(cds_dan, color_cells_by = "pseudotime",
           label_branch_points = FALSE, label_roots = FALSE,
           label_leaves = FALSE, cell_size = 1)

# Density plot across disease conditions
pseud_dan <- pseudotime(cds_dan)
pseud_dan <- pseud_dan[is.finite(pseud_dan)]
meta_dan_filt <- meta_dan[rownames(meta_dan) %in% names(pseud_dan), ]
meta_dan_filt$pseudotime <- pseud_dan

p_dens_dan <- ggplot(meta_dan_filt, aes(x = pseudotime, fill = Disease)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  labs(title = "DaN pseudotime by disease")
save_plot(p_dens_dan, here("figures/Figure_4/DaN_pseudotime_density.pdf"))

# switchDE along pseudotime
tpm_dan <- tpm_log(as.matrix(GetAssayData(sn_dan, slot = "counts")), gene_lengths)
tpm_dan_filt <- filter_expr(tpm_dan[, names(pseud_dan)])
sde_dan <- switchde(tpm_dan_filt, pseud_dan)
sde_dan <- dplyr::arrange(sde_dan, qval)
sde_dan <- sde_dan[sde_dan$qval < 0.01 & abs(sde_dan$k) > 0.03, ]
rm(tpm_dan, tpm_dan_filt); gc()

# ── Slingshot DaN — correlation with monocle3 pseudotime ─────────────────────
sce_dan <- as.SingleCellExperiment(sn_dan, assay = "RNA")

set.seed(SEED)
sce_dan <- slingshot(sce_dan, reducedDim = "UMAP", clusterLabels = "Disease")

meta_dan_filt$slingshot_ps <- sce_dan$slingPseudotime_1[
  match(rownames(meta_dan_filt), colnames(sce_dan))]

ps_corr <- cor.test(sce_dan$slingPseudotime_1, pseud_dan[colnames(sce_dan)],
                    use = "complete.obs")
message(sprintf("Monocle3 ~ slingshot correlation: r = %.3f, p = %.2e",
                ps_corr$estimate, ps_corr$p.value))

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: DaN_0 + DaN_1 — tradeSeq along slingshot trajectory
# ══════════════════════════════════════════════════════════════════════════════

sn_dan01 <- subset(sn_combined, subset = CellSubType %in% c("DaN_0", "DaN_1"))
sce_dan01 <- as.SingleCellExperiment(sn_dan01, assay = "RNA")

set.seed(SEED)
sce_dan01 <- slingshot(sce_dan01, reducedDim = "UMAP", clusterLabels = "Disease")

meta_dan01 <- data.frame(
  pseudotime = sce_dan01$slingPseudotime_1,
  disease    = sce_dan01$Disease,
  cell_type  = sce_dan01$CellSubType
)

genes_dan01 <- as.matrix(GetAssayData(sn_dan01, slot = "counts"))
set.seed(SEED)
sce_dan01_gam <- fitGAM(genes_dan01, sds = SlingshotDataSet(sce_dan01))
dan01_assoc   <- associationTest(sce_dan01_gam)
dan01_startend <- startVsEndTest(sce_dan01_gam)

# Heatmap: genes changing from start to end (|logFC| > 1, p < 0.05)
ps_intervals <- cut(sce_dan01$slingPseudotime_1, breaks = 3,
                    labels = c("ps_int1", "ps_int2", "ps_int3"))
names(ps_intervals) <- colnames(sn_dan01)
sn_dan01$ps_interval <- ps_intervals

sig_down <- rownames(dan01_startend)[
  dan01_startend$pvalue < 0.05 & dan01_startend$logFClineage1 < -1]
sig_up <- rownames(dan01_startend)[
  dan01_startend$pvalue < 0.05 & dan01_startend$logFClineage1 > 1]

avg_exp_ps <- AverageExpression(sn_dan01, group.by = "ps_interval",
                                features = c(sig_up, sig_down))
avg_exp_mat <- rbind(avg_exp_ps$RNA[sig_up, ], avg_exp_ps$RNA[sig_down, ])

pheatmap::pheatmap(log(avg_exp_mat + 1),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   scale = "row", fontsize = 2,
                   filename = here("figures/Figure_4/DaN01_pseudotime_heatmap.pdf"))

rm(genes_dan01); gc()

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: DaN_1 → DaN_3 cell state transition (CTR only, switchDE)
# ══════════════════════════════════════════════════════════════════════════════

sn_dan13_ctr <- subset(sn_combined,
                        subset = CellSubType %in% c("DaN_1", "DaN_3") &
                          Disease == CTRL_LABEL)
sce_dan13 <- as.SingleCellExperiment(sn_dan13_ctr, assay = "RNA")

set.seed(SEED)
sce_dan13 <- slingshot(sce_dan13, reducedDim = "UMAP",
                       clusterLabels = "CellSubType")

tpm_dan13  <- tpm_log(as.matrix(GetAssayData(sn_dan13_ctr, slot = "counts")),
                      gene_lengths)
tpm_dan13f <- filter_expr(tpm_dan13)

sde_dan13 <- switchde(tpm_dan13f, sce_dan13$slingPseudotime_1)
sde_dan13 <- dplyr::arrange(sde_dan13, qval)
sde_dan13 <- sde_dan13[sde_dan13$qval < 0.05, ]
dan_1and3_ctr_switchde <- as.data.frame(sde_dan13)
rm(tpm_dan13, tpm_dan13f); gc()

# PPI gene modules along three pseudotime intervals
ppi_net <- read.delim(here("data/external/ppi_net_2023.tsv"), header = TRUE)

split_and_build_ppi <- function(switchde_df, t_lo, t_hi, ppi_net) {
  sub_genes <- switchde_df$gene[
    (t_lo == -Inf | switchde_df$t0_sc > t_lo) &
      (t_hi == Inf | switchde_df$t0_sc <= t_hi)]

  subnet <- ppi_net[ppi_net[, 1] %in% sub_genes & ppi_net[, 2] %in% sub_genes, ]
  if (nrow(subnet) == 0) return(NULL)

  vertices <- unique(c(as.character(subnet[, 1]), as.character(subnet[, 2])))
  net  <- graph_from_data_frame(d = subnet, vertices = vertices, directed = FALSE)
  cl   <- cluster_louvain(net)
  data.frame(id = names(membership(cl)), module = as.integer(membership(cl)))
}

ppi_modules_t1 <- split_and_build_ppi(dan_1and3_ctr_switchde, -Inf,  9,   ppi_net)
ppi_modules_t2 <- split_and_build_ppi(dan_1and3_ctr_switchde,  9,   11,   ppi_net)
ppi_modules_t3 <- split_and_build_ppi(dan_1and3_ctr_switchde, 11,   Inf,  ppi_net)

# Average expression of RAS pathway genes across pseudotime intervals
ras <- read.table(here("data/external/ras_genes.txt"), header = FALSE)

sn_dan13_ctr$ps_range <- cut(sce_dan13$slingPseudotime_1,
                              breaks = PS_BREAKS,
                              include.lowest = TRUE)

avg_ras <- AverageExpression(sn_dan13_ctr,
                             features = ras[, 1],
                             group.by = "ps_range")
avg_ras_melt <- reshape2::melt(as.data.frame(avg_ras$RNA))

p_ras <- ggplot(avg_ras_melt, aes(x = variable, y = Var1, fill = value)) +
  geom_tile(color = "white", lwd = 1.2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = median(avg_ras_melt$value, na.rm = TRUE)) +
  labs(x = "Pseudotime interval", y = "Gene", fill = "Avg. expr.") +
  theme_minimal()
save_plot(p_ras, here("figures/Figure_4/DaN13_RAS_pseudotime.pdf"), width = 8, height = 6)

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: ODC pseudotime (using already-integrated data)
# ══════════════════════════════════════════════════════════════════════════════

# Subset from the integrated atlas (avoids re-processing and batch effects)
sn_odc <- subset(sn_combined, subset = CellType == "ODC")
message("ODC cells (from integrated atlas): ", ncol(sn_odc))

set.seed(SEED)
sn_odc <- NormalizeData(sn_odc, verbose = FALSE)
sn_odc <- FindVariableFeatures(sn_odc, verbose = FALSE)
sn_odc <- ScaleData(sn_odc, verbose = FALSE)
sn_odc <- RunPCA(sn_odc, verbose = FALSE)
sn_odc <- FindNeighbors(sn_odc, dims = 1:30, verbose = FALSE)

set.seed(SEED)
sn_odc <- FindClusters(sn_odc, resolution = 0.6, verbose = FALSE)
sn_odc <- RunUMAP(sn_odc, dims = 1:30, seed.use = SEED, verbose = FALSE)

sce_odc <- as.SingleCellExperiment(sn_odc, assay = "RNA")

set.seed(SEED)
sce_odc <- slingshot(sce_odc, reducedDim = "UMAP", clusterLabels = "Disease")

meta_odc2 <- data.frame(
  pseudotime = sce_odc$slingPseudotime_1,
  disease    = sce_odc$Disease,
  cell_type  = sce_odc$CellSubType
)
meta_odc2_filt <- meta_odc2[meta_odc2$cell_type == "ODC_2", ]

p_odc_dens <- ggplot(meta_odc2_filt, aes(x = pseudotime, fill = disease)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  labs(title = "ODC_2 pseudotime by disease")
save_plot(p_odc_dens, here("figures/Figure_5/ODC2_pseudotime_density.pdf"))

genes_odc  <- as.matrix(GetAssayData(sn_odc, slot = "counts"))
set.seed(SEED)
sce_odc_gam      <- fitGAM(genes_odc, sds = SlingshotDataSet(sce_odc))
odc2_assoc       <- associationTest(sce_odc_gam)
odc2_startend    <- startVsEndTest(sce_odc_gam)
rm(genes_odc); gc()

# ── Save ──────────────────────────────────────────────────────────────────────
ensure_dirs(here("results/pseudotime"))

save(cds_dan, sde_dan, sce_dan, ps_corr,
     dan01_assoc, dan01_startend, sce_dan01_gam,
     dan_1and3_ctr_switchde, ppi_modules_t1, ppi_modules_t2, ppi_modules_t3,
     file = here("results/pseudotime/ps_dan.RData"))

save(sce_odc, sce_odc_gam, odc2_startend,
     file = here("results/pseudotime/ps_odc.RData"))

message("All pseudotime results saved.")
log_session("07_pseudotime_analysis")
