# 02_cell_type_annotation.R
# Manual cell-type annotation using canonical marker genes.
# Project: deep_snRNAseqAtlas_humanSN — Volpato et al. (2026)

library(Seurat)
library(here)
source(here("scripts/utils.R"))

# ── Load integrated data ──────────────────────────────────────────────────────
load(here("data/processed/sn_atlas_integrated.RData"))
message("Loaded: ", ncol(sn_combined), " cells, ", nlevels(Idents(sn_combined)), " clusters")

# ── Marker genes for dot-plot QC ─────────────────────────────────────────────
markers_all <- c(
  "GFRA2", "CALCR", "CRYM", "CCDC68", "PPP1R17", "TRHR",
  "TH", "GRIN2C", "SLC18A2", "DCX", "LMX1A",
  "GAD2", "GAD1", "NTSR1", "SLC6A3", "KCNJ6", "PITX3",
  "CHRNA4", "LMX1B", "GRIK3", "ALDH1A1",
  "AQP4", "GFAP", "VCAN", "MOBP", "MOG",
  "CSF1R", "CD8A", "PTPRC", "CLDN5", "PTH1R"
)

DefaultAssay(sn_combined) <- "RNA"

p_dot <- DotPlot(sn_combined, features = markers_all, dot.scale = 8) + RotatedAxis()
save_plot(p_dot, here("figures/Figure_1/dotplot_markers_all.pdf"), width = 14, height = 5)

# ── Cluster → cell-type mapping ───────────────────────────────────────────────
# NOTE: update this vector if cluster numbering changes after re-running script 01
cell_type_map <- c(
  "ODC", "ODC", "DaN", "DaN", "Microglia",
  "Astrocyte", "OPC", "GABA", "DaN", "DaN", "DaN", "Tcell"
)
names(cell_type_map) <- levels(sn_combined)

# Validate mapping length matches number of clusters
stopifnot(
  "cell_type_map length must equal number of clusters" =
    length(cell_type_map) == nlevels(sn_combined)
)

sn_combined$cluster  <- Idents(sn_combined)
sn_combined          <- RenameIdents(sn_combined, cell_type_map)
sn_combined$CellType <- Idents(sn_combined)

message("Cell types assigned: ", paste(sort(unique(sn_combined$CellType)), collapse = ", "))
message("Cells per type:")
print(table(sn_combined$CellType))

# ── UMAP split by disease ─────────────────────────────────────────────────────
p_umap <- DimPlot(sn_combined, reduction = "umap",
                  split.by = "Disease", group.by = "CellType")
save_plot(p_umap, here("figures/Figure_1/umap_celltypes_by_disease.pdf"), width = 16, height = 5)

# ── Save ──────────────────────────────────────────────────────────────────────
save(sn_combined, file = here("data/processed/sn_atlas_annotated.RData"))
message("Saved annotated object.")

log_session("02_cell_type_annotation")
