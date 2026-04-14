```{r scenic, echo=TRUE}
### Initialize settings

suppressPackageStartupMessages({
   library(SCENIC)
   library(AUCell)
   library(RcisTarget)
   library(SCopeLoomR)
   library(KernSmooth)
   library(BiocParallel)
   library(ggplot2)
   library(data.table)
   library(grid)
   library(ComplexHeatmap)
   library(Seurat)
   library(doRNG)
 })


load(here("data/processed/sn_atlas_annotated_subtype.RData"))
sn_combined_tmp <- subset(sn_combined,subset=CellType=="DaN")
sn_combined_tmp <- subset(sn_combined_tmp,subset=Disease=="CTR")


genes_dan_ctr<-as.matrix(GetAssayData(sn_combined_tmp, slot = "counts"))
metadata_danALLctr<-as.data.frame(tmp@meta.data)
org <- "hgnc"
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
dbDir <-"/scratch/c.mpmvv/Postmortem_02_NovaSeq/src/"

scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle="dans", nCores=10)
scenicOptions@inputDatasetInfo$cellInfo<-metadata_danALLctr
scenicOptions@inputDatasetInfo$colVars<-metadata_danALLctr$CellSubType

data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9

### Co-expression network

exprMat_filtered <- genes_dan_ctr
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions_1.Rds")
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions_2.Rds")
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log)
saveRDS(scenicOptions, file="int/scenicOptions_3.Rds")


Idents(sn_combined_tmp)<-sn_atlas_dans_red$CellSubType
cellInfo <- data.frame(seuratCluster=Idents(sn_combined_tmp))
rownames(cellInfo)<-colnames(sn_combined_tmp)

regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo[,1]),
                                     function(cells) rowMeans(getAUC(regulonAUC_dans)[,cells]))




regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",row_names_gp =gpar(fontsize = 5))

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]

                                     
