# 10_coexpression_network.R
# BigScale2 co-expression network for DaN_3, with Louvain community detection.
# Project: deep_snRNAseqAtlas_humanSN — Volpato et al. (2026)

library(Seurat)
library(BigScale2)
library(igraph)
library(here)
source(here("scripts/utils.R"))

# ── Parameters ────────────────────────────────────────────────────────────────
CORR_THRESHOLD <- 0.8   # |r| cut-off for keeping co-expression edges
TARGET_SUBTYPE <- "DaN_3"

# ── Load data ─────────────────────────────────────────────────────────────────
load(here("data/processed/sn_atlas_annotated_subtype.RData"))

sn_dan3 <- subset(sn_combined, subset = CellSubType == TARGET_SUBTYPE)
message("Cells in ", TARGET_SUBTYPE, ": ", ncol(sn_dan3))

# ── Co-expression network (BigScale2) ─────────────────────────────────────────
# Note: BigScale2 works on raw count data
data_dan3 <- as.matrix(GetAssayData(sn_dan3, slot = "counts"))

results_dan3 <- compute.network(
  expr.data  = data_dan3,
  gene.names = rownames(data_dan3),
  clustering = "direct"
)
rm(data_dan3); gc()

# ── Extract upper-triangle edges ──────────────────────────────────────────────
corr_mat <- as.data.frame(results_dan3$correlations)

coexp_dan3 <- data.frame(
  gene_a = rownames(corr_mat)[row(corr_mat)[upper.tri(corr_mat)]],
  gene_b = colnames(corr_mat)[col(corr_mat)[upper.tri(corr_mat)]],
  corr   = corr_mat[upper.tri(corr_mat)]
)
rm(corr_mat); gc()

# Apply correlation threshold
coexp_dan3_filt <- coexp_dan3[abs(coexp_dan3$corr) > CORR_THRESHOLD, ]
message(sprintf("Edges with |r| > %.1f: %d", CORR_THRESHOLD, nrow(coexp_dan3_filt)))

# ── Louvain community detection ───────────────────────────────────────────────
coexp_genes <- unique(c(coexp_dan3_filt$gene_a, coexp_dan3_filt$gene_b))

net_dan3  <- graph_from_data_frame(d = coexp_dan3_filt,
                                   vertices = coexp_genes,
                                   directed = FALSE)
cl_net    <- cluster_louvain(net_dan3)
memb_dan3 <- as.data.frame(membership(cl_net))

memb_dan3_df <- data.frame(
  id     = rownames(memb_dan3),
  module = memb_dan3[, 1]
)
message("Co-expression modules found: ", length(unique(memb_dan3_df$module)))

# ── Save ──────────────────────────────────────────────────────────────────────
ensure_dirs(here("results/coexpression"))
save(results_dan3, coexp_dan3_filt, memb_dan3_df,
     file = here("results/coexpression/coexp_DaN3.RData"))
message("Co-expression network results saved.")

log_session("10_coexpression_network")
