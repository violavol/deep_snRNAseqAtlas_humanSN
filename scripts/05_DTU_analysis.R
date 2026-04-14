### Isoform Analysis 

library(slingshot)
library(fishpond)
library(scran)
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://useast.ensembl.org")
library(tximeta)

coldata<-read.delim("data/processed/coldata_ALL_CTRandPD56_modified_nfs",h=T,stringsAsFactors=F)

suppressPackageStartupMessages(library(SummarizedExperiment))
y <- tximeta(coldata,dropInfReps=T)

# diff transcript express

y <- labelKeep(y, minCount = 3, minN = 10)
y <- y[mcols(y)$keep,]
set.seed(1)
assays(y) <- lapply(assays(y), as.matrix) # make dense matrices
y <- scaleInfReps(y, lengthCorrect=FALSE, sfFun= computeSumFactors)
y <- swish(y, x="condition",  quiet=TRUE)

# diff isoform usage

iso <- isoformProportions(y)
iso <- swish(iso, x="condition",nperms=64)

df<-mcols(iso)[,c("log2FC","qvalue","gene","tx_id")]

load("data/processed/Transcript_expression_allCelltype_CTRandPD56.RData") # to get all data and get number of DTUs

# OPCs to ODCs differentiation using isoform data

load("data/processed/pseudotime_isoforms_OPC_ODCs_CTRandPD56.RData")

isoform_OPC_ODC2_CTRandPD56<-as.matrix(assays(y)[[2]])
meta_isoform_OPC_ODC2_CTRandPD56<-data.frame(celltype=y$cellSubType,disease=y$condition,cell=y$Barcode,sample=y$names)
meta_isoform_OPC_ODC2_CTRandPD56$names<-paste(meta_isoform_OPC_ODC2_CTRandPD56$cell,meta_isoform_OPC_ODC2_CTRandPD56$sample,sep="_")
colnames(isoform_OPC_ODC2_CTRandPD56)<-meta_isoform_OPC_ODC2_CTRandPD56$names
rownames(meta_isoform_OPC_ODC2_CTRandPD56)<-meta_isoform_OPC_ODC2_CTRandPD56$names

chip_iso_opc_odc2_ctrpd <- CreateSeuratObject(counts = isoform_OPC_ODC2_CTRandPD56, project = "isoform", min.cells = 0, min.features = 0,meta.data= meta_isoform_OPC_ODC2_CTRandPD56)

# done on controls and pd samples separately
chip_iso_opc_odc2_ctr <- subset(chip_iso_opc_odc2_ctrpd,subset=disease=="CTR")
chip_iso_opc_odc2_ctr <- NormalizeData(chip_iso_opc_odc2_ctr)
chip_iso_opc_odc2_ctr <- FindVariableFeatures(chip_iso_opc_odc2_ctr, selection.method = "vst", nfeatures = 2000)

chip_iso_opc_odc2_ctr <- ScaleData(chip_iso_opc_odc2_ctr)
chip_iso_opc_odc2_ctr <- RunPCA(chip_iso_opc_odc2_ctr, features = VariableFeatures(object = chip_iso_opc_odc2_ctrpd))
chip_iso_opc_odc2_ctr <- FindNeighbors(chip_iso_opc_odc2_ctr, dims = 1:10)
chip_iso_opc_odc2_ctr <- FindClusters(chip_iso_opc_odc2_ctr, resolution = 0.6)
chip_iso_opc_odc2_ctr <- RunUMAP(chip_iso_opc_odc2_ctr, dims = 1:10)

sce_iso_opc_odc2_ctr <- as.SingleCellExperiment(chip_iso_opc_odc2_ctr, assay = "RNA")
sce_iso_opc_odc2_ctr <- slingshot(sce_iso_opc_odc2_ctr, reducedDim = 'UMAP', clusterLabels = 'seurat_clusters')
meta_isoform_OPC_ODC2_CTR<-meta_isoform_OPC_ODC2_CTRandPD56[meta_isoform_OPC_ODC2_CTRandPD56$disease=="CTR",]
meta_isoform_OPC_ODC2_CTR$ps<-sce_iso_opc_odc2_ctr$slingPseudotime_1

ggplot(meta_isoform_OPC_ODC2_CTR, aes(x=ps, fill=celltype)) +
    geom_density(alpha=.5) + theme_bw()

save(meta_isoform_OPC_ODC2_CTR,meta_isoform_OPC_ODC2_PD_B56, file= here("results/pseudotime/OPC_ODCs_DTU_ps.RData"))


isoform_OPC_ODC2_CTR<-as.matrix(GetAssayData(chip_iso_opc_odc2_ctr, slot = "counts"))
isoform_OPC_ODC2_CTR<-log(isoform_OPC_ODC2_CTR+1)
sde_iso_opc_odc2_ctr <- switchde(isoform_OPC_ODC2_CTR, sce_iso_opc_odc2_ctr$slingPseudotime_1)

# test enrichment of pd risk isoforms in OPCs and ODC_2 along differentiation trajectory:
pd_risk_opc_isoform<-read.table("pd_50k_gene_and_iso_opc.gsa.sets.genes.out",h=T)
pd_risk_odc2_isoform<-read.table("pd_50k_gene2000_and_iso_odc2_only.gsa.sets.genes.out",h=T)




