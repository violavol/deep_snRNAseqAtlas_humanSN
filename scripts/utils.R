# utils.R — Shared utility functions
# Project: deep_snRNAseqAtlas_humanSN
# Volpato et al. (2026)

library(ggplot2)
library(dplyr)
library(Seurat)

# ── Jaccard index between two gene sets ───────────────────────────────────────
jaccard <- function(a, b) {
  inter <- length(intersect(a, b))
  union <- length(a) + length(b) - inter
  inter / union
}

# ── Save a ggplot to file ─────────────────────────────────────────────────────
save_plot <- function(plot_obj, filename, width = 6, height = 4) {
  ggsave(filename, plot = plot_obj, width = width, height = height)
}

# ── DotPlot wrapper with rotated axis ─────────────────────────────────────────
dotplot_rot <- function(seurat_obj, features, cluster_ids = NULL, scale = 8) {
  DotPlot(seurat_obj, features = features, idents = cluster_ids, dot.scale = scale) +
    RotatedAxis() +
    theme_minimal()
}

# ── W-test (signed-rank linear model) for cell-type specificity ───────────────
# Computes a specificity statistic for each gene × cell-type pair using a
# linear model on signed ranks.  Parallelised over genes with BiocParallel.
#
# Args:
#   seurat_obj : Seurat object (raw counts used)
#   cell_col   : metadata column holding cell-type labels
#   subset_col : optional metadata column to pre-filter rows
#   subset_val : value(s) to keep in subset_col
#
# Returns:
#   data.frame (genes × cell-types) of t-statistics
run_wtest <- function(seurat_obj,
                      cell_col,
                      subset_col = NULL,
                      subset_val = NULL) {

  library(BiocParallel)

  signed_rank <- function(x) sign(x) * rank(abs(x))

  meta <- as.data.frame(seurat_obj@meta.data)

  # Optional pre-filter
  if (!is.null(subset_col) && !is.null(subset_val)) {
    meta <- meta[meta[[subset_col]] %in% subset_val, , drop = FALSE]
  }

  # Keep counts sparse; only materialise columns we need
  counts <- GetAssayData(seurat_obj, slot = "counts")[, rownames(meta), drop = FALSE]
  cell_types <- sort(unique(meta[[cell_col]]))
  genes      <- rownames(counts)

  message(sprintf(
    "run_wtest: %d genes × %d cell types (%d cells)",
    length(genes), length(cell_types), ncol(counts)
  ))

  compute_one_gene <- function(gene) {
    y <- as.numeric(counts[gene, ])
    yr <- signed_rank(y)
    vapply(cell_types, function(ct) {
      x <- ifelse(meta[[cell_col]] == ct, 1L, -1L)
      df <- data.frame(yr = yr, x = x)
      coef(summary(lm(yr ~ x, data = df)))[2, 3]
    }, numeric(1))
  }

  results <- bplapply(genes, compute_one_gene,
                      BPPARAM = MulticoreParam(workers = max(1L, parallel::detectCores() - 1L)))

  mat <- do.call(rbind, results)
  rownames(mat) <- genes
  colnames(mat) <- cell_types
  as.data.frame(mat)
}

# ── Log session info to file ──────────────────────────────────────────────────
log_session <- function(script_name, log_dir = here::here("logs")) {
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  out_path <- file.path(log_dir, paste0("sessionInfo_", script_name, ".txt"))
  writeLines(capture.output(sessionInfo()), out_path)
  message("Session info written to: ", out_path)
}

# ── Helper: ensure output directories exist ───────────────────────────────────
ensure_dirs <- function(...) {
  dirs <- c(...)
  invisible(lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE))
}
