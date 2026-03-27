library(Seurat)
library(here)

# Load annotated object
load(here("data/processed/sn_atlas_annotated.RData"))

# Sub-cluster dopaminergic neurons
DA_cells <- subset(sn_combined, idents = "DaN")
meta<-as.data.frame(DA_cells@meta.data)
meta<-meta[meta$Sample_v2!="14_133",] # outlier sample in both UMAP plot of all DaNs and in pseudotime analysis 
DA_cells <- DA_cells[,colnames(DA_cells)%in%rownames(meta)]

tmp_cell <- DA_cells
tmp_cell.list <- SplitObject(tmp_cell, split.by = "Disease")

tmp_cell.list <- lapply(X = tmp_cell.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = tmp_cell.list)
tmp_cell.list <- lapply(X = tmp_cell.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

tmp_cell.anchors <- FindIntegrationAnchors(object.list = tmp_cell.list, anchor.features = features, reduction = "rpca", k.anchor=40)
tmp_cell.combined <- IntegrateData(anchorset = tmp_cell.anchors) 

#tmp_cell.combined <- IntegrateData(anchorset = tmp_cell.anchors) # for odc
#tmp_cell.combined <- IntegrateData(anchorset = tmp_cell.anchors,k.weight=80) # for microglia
#tmp_cell.combined <- IntegrateData(anchorset = tmp_cell.anchors,k.weight=120) # for astrocytes

DefaultAssay(tmp_cell.combined) <- "integrated"

tmp_cell.combined <- ScaleData(tmp_cell.combined, verbose = FALSE)
tmp_cell.combined <- RunPCA(tmp_cell.combined, npcs = 40, verbose = FALSE)
tmp_cell.combined <- RunUMAP(tmp_cell.combined, reduction = "pca", dims = 1:40)
tmp_cell.combined <- FindNeighbors(tmp_cell.combined, reduction = "pca", dims = 1:40)
tmp_cell.combined <- FindClusters(tmp_cell.combined, resolution = 0.2)

DimPlot(tmp_cell.combined, reduction = "umap", split.by = "Disease",group.by="seurat_clusters")

DefaultAssay(tmp_cell.combined) <- "RNA"
dan_markers <- c("SOX6", "GRIA3", "DCX", "AGTR1", "LMX1B", "GRIK3", "RET", "GFRA2", "PITX3", "CHRNA4", "SLC18A2", "SLC6A3", "TH", "KCNJ6", "ALDH1A1", "TMEM255A", "LGI1",   "TMEFF2")

DotPlot(tmp_cell.combined, features = dan_markers, dot.scale = 8,idents = c(0,1,2,3)) + # clusters 4 and 5 are not used as too small and only present in ILBD_B3-4
    RotatedAxis()


# markers_ODC=c("PLXDC2","PLP1","SPARC","DHCR24","TUBA1A","PMP2","RBFOX1","AFF3","FMN1","PALM2","HHIP","OPALIN","LAMA2") # selected from https://www.biorxiv.org/content/10.1101/2022.03.22.485367v1.full.pdf

# markers_microglia=c("PTPRC","ITGAM","AIF1","C1QA","CTSS","CD14","CSF3R","ARGLU1","FAM46A","IFIT3","ISG15","MRC1","TNF","CD83","EGR2","TNFSF18","CCL8","TFRC","PCNA","RASGEF1C","AC008691.1", "TLN2","GPNMB","ACSL1","CXCR4","NTM", "MAGI2", "SCD", "PLP1","NRG3")  (check https://www.biorxiv.org/content/10.1101/2022.03.22.485367v1.full.pdf and https://www.nature.com/articles/s41467-020-19737-2) markers.to.plot=c("MRC1","IL10","ABCC4","CSF2RA","CSF3R","TFRC","KLF4","PTGS1","DOCK8","KCNQ1","PTPRC","GPNMB","SPP1","TYROBP","TREM2","TLR2","MS4A4A","IL13RA1","INPP5D","ITGAM","ADAP2","APBB1IP","SP140L","VAV1") # macrophage M2: "MRC1","IL10","ABCC4","CSF2RA","CSF3R","TFRC","KLF4","PTGS1","DOCK8"

# markers_astrocyte=c("GABRA2","EDNRB","PPFIA2","KCNJ16","BHLHE40","SLC6A11","PAK3","GRM5","PTPRT","FAT3","EPHA6","EPHB1","ADGRV1","SPOCK1","SPSB1","MYO1E","SLC24A2","ELMO1","NKAIN2","PLP1","S100B","TMSB4X") 

