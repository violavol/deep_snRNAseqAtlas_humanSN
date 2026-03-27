library(Seurat)
library(here)

# Load integrated data
load(here("data/processed/sn_atlas_integrated.RData"))

# Cell type annotation
DefaultAssay(sn_combined) <- "RNA"

markersALL<-c("GFRA2","CALCR","CRYM","CCDC68","PPP1R17","TRHR","TH","GRIN2C","SLC18A2","DCX","LMX1A","GAD2","GAD1","NTSR1","SLC6A3","KCNJ6","PITX3","CHRNA4","LMX1B","GRIK3","ALDH1A1","AQP4","GFAP","VCAN","MOBP","MOG","CSF1R","CD8A","PTPRC","CLDN5","PTH1R")

DotPlot(sn_combined, features = markersALL, dot.scale = 8) +
    RotatedAxis()

sn_combined$cluster <- Idents(sn_combined)
cellType <- c("ODC","ODC","DaN","DaN","Microglia","Astrocyte","OPC","GABA","DaN","DaN","DaN","Tcell")
names(cellType) <- levels(sn_combined)
sn_combined <- RenameIdents(sn_combined, cellType)
sn_combined$CellType<-Idents(sn_combined)

DimPlot(sn_combined, reduction = "umap", split.by = "Disease",group.by="CellType")

# Save annotated object
save(sn_combined, file = here("data/processed/sn_atlas_annotated.RData"))
