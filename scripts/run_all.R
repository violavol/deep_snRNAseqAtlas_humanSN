# run_all.R
# Master script — runs the full analysis pipeline in order.
# Project: deep_snRNAseqAtlas_humanSN — Volpato et al. (2026)
#
# Usage:
#   Rscript scripts/run_all.R
#   # or interactively, source() each step you need

library(here)

scripts <- c(
  "01_data_preprocessing.R",
  "02_cell_type_annotation.R",
  "03_subtype_analysis.R",
  "04_DEG_analysis.R",
  "05_DTU_analysis.R",
  "06_pathway_PPI_analysis.R",
  "07_pseudotime_analysis.R",
  "08_cell_proportions.R",
  "09_cellcommunication_analysis.R",
  "10_coexpression_network.R",
  "11_SCENIC.R"
)

log_dir <- here("logs")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

for (script in scripts) {
  path <- here("scripts", script)
  message("\n", strrep("=", 60))
  message("▶  Running: ", script)
  message(strrep("=", 60))
  tryCatch(
    source(path, echo = FALSE),
    error = function(e) {
      message("✖  ERROR in ", script, ":\n", conditionMessage(e))
      stop(e)
    }
  )
  message("✔  Done: ", script)
}

message("\n✅ Full pipeline complete.")
