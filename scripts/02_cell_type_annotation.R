library(Seurat)
library(here)

# Load integrated data
load(here("data/processed/sn_atlas_integrated.RData"))

# Cell type annotation
DefaultAssay(combined_all) <- "RNA"

markersALL<-c("GFRA2","CALCR","CRYM","CCDC68","PPP1R17","TRHR","TH","GRIN2C","SLC18A2","DCX","LMX1A","GAD2","GAD1","NTSR1","SLC6A3","KCNJ6","PITX3","CHRNA4","LMX1B","GRIK3","ALDH1A1","AQP4","GFAP","VCAN","MOBP","MOG","CSF1R","CD8A","PTPRC","CLDN5","PTH1R")


for (cell_type in names(markers)) {
  sn_combined <- AddModuleScore(sn_combined, features = list(markers[[cell_type]]), name = cell_type)
}

# Save annotated object
save(sn_combined, file = here("data/processed/sn_atlas_annotated.RData"))
