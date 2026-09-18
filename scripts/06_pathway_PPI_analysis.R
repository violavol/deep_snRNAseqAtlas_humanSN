# 06_pathway_PPI_analysis.R
# GO enrichment on DEGs and PPI network community detection for DaN_3, ODC_2, OPC.
# Project: deep_snRNAseqAtlas_humanSN — Volpato et al. (2026)

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(igraph)
library(here)
source(here("scripts/utils.R"))

# ── Parameters ────────────────────────────────────────────────────────────────
PADJ_CUTOFF  <- 0.05
PPI_TOP_N    <- 1000    # top W-test genes per cell type for PPI network
GO_P_CUTOFF  <- 0.05
GO_Q_CUTOFF  <- 0.10

# ── Load results ──────────────────────────────────────────────────────────────
load(here("results/DEG_tables/DEG_results.RData"))   # deg_results, df_wtest_level1, df_wtest_level2_DaN

# ── GO enrichment helper ──────────────────────────────────────────────────────
go_enrichment <- function(gene_list, ont = "BP") {
  if (length(gene_list) == 0) return(NULL)

  gene_ids <- bitr(gene_list, fromType = "SYMBOL",
                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (nrow(gene_ids) == 0) return(NULL)

  enrichGO(
    gene          = gene_ids$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = GO_P_CUTOFF,
    qvalueCutoff  = GO_Q_CUTOFF
  )
}

# ── GO on DEG results ─────────────────────────────────────────────────────────
go_results_degs <- lapply(names(deg_results), function(st) {
  sig_genes <- rownames(deg_results[[st]])[deg_results[[st]]$p_val_adj < PADJ_CUTOFF]
  message(st, ": ", length(sig_genes), " significant DEGs")
  go_enrichment(sig_genes)
})
names(go_results_degs) <- names(deg_results)

# ── PPI network and community detection ──────────────────────────────────────
ppi_net <- read.delim(here("data/external/ppi_net_2023.tsv"), header = FALSE)
names(ppi_net)[1:2] <- c("gene_a", "gene_b")

build_ppi_subnetwork <- function(wtest_df, cell_type, top_n, ppi_net) {
  top_genes <- rownames(wtest_df)[order(-wtest_df[[cell_type]])[seq_len(top_n)]]

  subnet <- ppi_net[
    ppi_net$gene_a %in% top_genes &
      ppi_net$gene_b %in% top_genes, ]

  if (nrow(subnet) == 0) {
    warning("No PPI edges for ", cell_type, " — returning NULL")
    return(NULL)
  }

  vertices <- unique(c(subnet$gene_a, subnet$gene_b))
  net      <- graph_from_data_frame(d = subnet, vertices = vertices, directed = FALSE)
  cl       <- cluster_louvain(net)

  data.frame(
    id     = names(membership(cl)),
    module = as.integer(membership(cl)),
    cell_type = cell_type
  )
}

# Run for the three cell types with PD genetic risk convergence
target_cell_types <- c("DaN_3", "ODC_2", "OPC")
wtest_combined <- cbind(df_wtest_level1, df_wtest_level2_DaN)

ppi_membership <- lapply(target_cell_types, function(ct) {
  if (!ct %in% colnames(wtest_combined)) {
    warning("W-test column not found for ", ct)
    return(NULL)
  }
  message("Building PPI network for: ", ct)
  build_ppi_subnetwork(wtest_combined, ct, PPI_TOP_N, ppi_net)
})
names(ppi_membership) <- target_cell_types

# ── GO on PPI modules ─────────────────────────────────────────────────────────
go_results_ppi <- lapply(names(ppi_membership), function(ct) {
  memb <- ppi_membership[[ct]]
  if (is.null(memb)) return(NULL)

  lapply(sort(unique(memb$module)), function(mod) {
    genes <- memb$id[memb$module == mod]
    go_enrichment(genes)
  })
})
names(go_results_ppi) <- target_cell_types

# ── Save ──────────────────────────────────────────────────────────────────────
ensure_dirs(here("results/GO_PPI"))
save(go_results_degs, ppi_membership, go_results_ppi,
     file = here("results/GO_PPI/GO_PPI_results.RData"))
message("GO and PPI results saved.")

log_session("06_pathway_PPI_analysis")
