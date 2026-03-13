library(Seurat)
library(here)

# Load integrated data
load(here("data/processed/sn_atlas_integrated.RData"))

# Cell type annotation
# Example: based on canonical markers
markers <- list(
  Dopaminergic = c("TH","SLC6A3"),
  Astrocytes = c("GFAP","AQP4"),
  Microglia = c("CX3CR1","P2RY12")
)

for (cell_type in names(markers)) {
  sn_combined <- AddModuleScore(sn_combined, features = list(markers[[cell_type]]), name = cell_type)
}

# Save annotated object
save(sn_combined, file = here("data/processed/sn_atlas_annotated.RData"))
