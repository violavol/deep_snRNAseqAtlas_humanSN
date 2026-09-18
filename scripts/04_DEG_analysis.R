# 04_DEG_analysis.R
# Differential gene expression (CTR vs PD_B5-6) per cell sub-type, plus
# W-test specificity scores for MAGMA gene-set analysis.
# Project: deep_snRNAseqAtlas_humanSN — Volpato et al. (2026)

library(Seurat)
library(dplyr)
library(here)
source(here("scripts/utils.R"))

# ── Parameters ────────────────────────────────────────────────────────────────
CTRL_LABEL    <- "CTR"
PD_LABEL      <- "PD_B5-6"
MAGMA_TOP_N   <- 500    # top N genes per cell type for MAGMA input
GENE_LOC_FILE <- here("data/external/NCBI37.3.gene.loc")

# ── Load data ─────────────────────────────────────────────────────────────────
load(here("data/processed/sn_atlas_annotated_subtype.RData"))

# ── 1. DEG analysis: CTR vs PD_B5-6 per sub-type ────────────────────────────
sn_deg <- subset(sn_combined, subset = Disease %in% c(CTRL_LABEL, PD_LABEL))
sn_deg$clust_disease <- paste(sn_deg$CellSubType, sn_deg$Disease, sep = "_")
Idents(sn_deg) <- "clust_disease"

sub_types <- sort(unique(sn_combined$CellSubType))
sub_types <- sub_types[!is.na(sub_types)]

deg_results <- lapply(sub_types, function(st) {
  ident_1 <- paste0(st, "_", PD_LABEL)
  ident_2 <- paste0(st, "_", CTRL_LABEL)

  # Guard: skip if either group is absent
  if (!ident_1 %in% levels(Idents(sn_deg)) ||
      !ident_2 %in% levels(Idents(sn_deg))) {
    message("Skipping ", st, " — one condition missing")
    return(NULL)
  }

  FindMarkers(sn_deg, ident.1 = ident_1, ident.2 = ident_2,
              verbose = FALSE)
})
names(deg_results) <- sub_types
deg_results <- Filter(Negate(is.null), deg_results)

message("DEG analysis complete for: ", paste(names(deg_results), collapse = ", "))

# ── 2. W-test: level-1 (all cell types) ──────────────────────────────────────
# Runs in parallel via BiocParallel (see run_wtest in utils.R)
message("Running W-test — level 1 (all cell types, CTR cells only)...")
df_wtest_level1 <- run_wtest(
  seurat_obj = sn_combined,
  cell_col   = "CellType",
  subset_col = "Disease",
  subset_val = CTRL_LABEL
)
write.table(df_wtest_level1, here("results/DEG_tables/Wtest_level1.tsv"),
            quote = FALSE, sep = "\t")

# ── 3. W-test: level-2 (DaN sub-types, CTR only) ─────────────────────────────
message("Running W-test — level 2 (DaN sub-types, CTR cells only)...")
dan_ctr <- subset(sn_combined,
                  subset = CellType == "DaN" & Disease == CTRL_LABEL)

df_wtest_level2_DaN <- run_wtest(
  seurat_obj = dan_ctr,
  cell_col   = "CellSubType"
)
write.table(df_wtest_level2_DaN, here("results/DEG_tables/Wtest_DaN_subtypes.tsv"),
            quote = FALSE, sep = "\t")

rm(dan_ctr); gc()

# ── 4. MAGMA input preparation ────────────────────────────────────────────────
gene_loc <- read.delim(GENE_LOC_FILE, header = FALSE)
names(gene_loc)[6] <- "id"

prepare_magma_input <- function(wtest_df, top_n, gene_loc, out_file) {
  cell_types <- colnames(wtest_df)

  gene_sets <- lapply(cell_types, function(ct) {
    top_genes <- rownames(wtest_df)[order(-wtest_df[[ct]])[seq_len(top_n)]]
    data.frame(cellt = ct, id = top_genes)
  })
  gene_sets_df <- do.call(rbind, gene_sets)

  merged <- merge(gene_sets_df, gene_loc[, c("id", "V1")], by = "id")

  lines <- vapply(cell_types, function(ct) {
    ids <- merged$V1[merged$cellt == ct]
    paste(c(ct, ids), collapse = " ")
  }, character(1))

  write.table(data.frame(lines), file = out_file,
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  message("MAGMA input written to: ", out_file)
}

prepare_magma_input(
  wtest_df  = df_wtest_level1,
  top_n     = MAGMA_TOP_N,
  gene_loc  = gene_loc,
  out_file  = here("results/DEG_tables/input_magma_level1.tsv")
)

prepare_magma_input(
  wtest_df  = df_wtest_level2_DaN,
  top_n     = MAGMA_TOP_N,
  gene_loc  = gene_loc,
  out_file  = here("results/DEG_tables/input_magma_DaN_subtypes.tsv")
)

# ── Save all results ──────────────────────────────────────────────────────────
ensure_dirs(here("results/DEG_tables"))
save(deg_results, df_wtest_level1, df_wtest_level2_DaN,
     file = here("results/DEG_tables/DEG_results.RData"))
message("All DEG results saved.")

log_session("04_DEG_analysis")
