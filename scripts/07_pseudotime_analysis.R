library(Seurat)
library(monocle3)
library(slingshot)
library(here)

# Pseudotime analysis for each cell type 

# DaNs:

load(here("data/processed/sn_atlas_annotated_subtype.RData"))

sn_combined_dan<-subset(sn_combined,subset=CellType=="DaN)
genesALL_DaN<-as.matrix(GetAssayData(sn_combined_dan, slot = "counts"))
metadataALL_DaN<-sn_combined_dan@meta.data

gene_annotation<-data.frame(gene_short_name=rownames(genesALL_DaN),genes=rownames(genesALL_DaN))
rownames(gene_annotation)<-rownames(genesALL_DaN)

cds_dan <- new_cell_data_set(genesALL_DaN,
                             cell_metadata = metadataALL_DaN,
                             gene_metadata = gene_annotation)

cds_dan <- preprocess_cds(cds_dan, num_dim = 3)

cds_dan <- reduce_dimension(cds_dan)
cds_dan <- cluster_cells(cds_dan, resolution=1e-5)
cds_dan <- learn_graph(cds_dan)

plot_cells(cds_dan, color_cells_by = "Disease",label_branch_points = F,label_roots = F,label_leaves = F, cell_size=1)

cds_dan <- order_cells(cds_dan)
plot_cells(cds_dan, color_cells_by = "pseudotime",label_branch_points = F,label_roots = F,label_leaves = F, cell_size=1)


# density plot

pseud_dan<-pseudotime(cds_dan)
pseud_dan<-pseud_dan[is.finite(pseud_dan)]
metadataALL_DaN_filt<-metadataALL_DaN[rownames(metadataALL_DaN)%in%names(pseud_dan),]
metadataALL_DaN_filt$pseudotime<-pseud_dan

ggplot(metadataALL_DaN_filt, aes(x=pseudotime, fill=Disease)) +
    geom_density(alpha=.5) + theme_bw()


# Find DE genes along pseudotime (switchDE)

prot_cod_len<-read.delim("data/processed/gene_length",h=T)
data<-as.matrix(GetAssayData(sn_combined_dan, slot = "counts"))
prot_cod_m<-prot_cod_len[match(rownames(data),prot_cod_len[,1]),]
data <- sweep(data, 1, STATS = prot_cod_m$Gene_Length/1000, FUN = "/")
data[is.na(data)]<-0
data <- sweep(data, 2, STATS = colSums(data)/(10^6), FUN = "/")
SN_TPM_log_red<-log(data+1)
X_filtered <- SN_TPM_log_red[rowMeans(SN_TPM_log_red) > 0.1 & rowMeans(SN_TPM_log_red > 0) > 0.2,]
sde_dan <- switchde(X_filtered, pseud_dan)
sde_dan <-arrange(sde_dan,qval)
sde_dan<-sde_dan[sde_dan$qval<0.01 & abs(sde_dan$k)>0.03,]


# Correlation with slingshot pseudotime:
# use sn_atlas_processed (with all cell types) and run Seurat classical pipeline, then subset only DaNs and run slingshot

sce_dans <- as.SingleCellExperiment(sn_combined_dan, assay = "RNA")
sce_dans <- slingshot(sce_dans, reducedDim = 'UMAP', clusterLabels = 'Disease')
meta_toplot<-data.frame(ps=sce_dans$slingPseudotime_1,disease=sce_dans$Disease,cellt=sce_dans$CellSubType)

ggplot(meta_toplot, aes(x=ps, fill=disease)) +
geom_density(alpha=.5) + theme_bw()

cor.test(sce_dans$slingPseudotime_1,pseudotime_DaNs_seuratIntegration_PC40res02[,2]) # corr = 0.61

# find DE genes along pseudotime (tradeseq)
# done for DaN_0 and DaN_1 subtypes across all disease condition

sn_combined_tmp<-subset(sn_combined,subset=CellSubType=="DaN_0" | CellSubType=="DaN_1")
sce_dans <- as.SingleCellExperiment(sn_combined_tmp, assay = "RNA")
sce_dans <- slingshot(sce_dans, reducedDim = 'UMAP', clusterLabels = 'Disease')
meta_toplot<-data.frame(ps=sce_dans$slingPseudotime_1,disease=sce_dans$Disease,cellt=sce_dans$CellSubType)

ggplot(meta_toplot, aes(x=ps, fill=disease)) +
     geom_density(alpha=.5) + theme_bw()

genes_dans<-as.matrix(GetAssayData(sn_combined_tmp, slot = "counts"))
sce_dans01_deg <- fitGAM(genes_dans,sds=SlingshotDataSet(sce_dans))
dans01_deg_tradeseq <- associationTest(sce_dans01_deg)
dans01_deg_tradeseqStartEnd <- startVsEndTest(sce_dans01_deg)

ps_dans01<-cut(sce_dans$slingPseudotime_1,breaks = 3,labels = c("ps_int1","ps_int2","ps_int3"))
names(ps_dans01)<-colnames(sn_atlas_processed_tmp)
sn_atlas_processed_tmp$ps_interval<-ps_dans01
avg_exp_ps_dan01_down<-AverageExpression(object = sn_atlas_processed_tmp, group.by = c('ps_interval'),features = rownames(dans01_deg_tradeseqStartEnd[dans01_deg_tradeseqStartEnd$pvalue<0.05 & dans01_deg_tradeseqStartEnd$logFClineage1< -1,]))
avg_exp_ps_dan01_up<-AverageExpression(object = sn_atlas_processed_tmp, group.by = c('ps_interval'),features = rownames(dans01_deg_tradeseqStartEnd[dans01_deg_tradeseqStartEnd$pvalue<0.05 & dans01_deg_tradeseqStartEnd$logFClineage1> 1,]))
avg_exp_ps_dan01_all<-rbind(avg_exp_ps_dan01_up$RNA,avg_exp_ps_dan01_down$RNA)

pheatmap::pheatmap(log(avg_exp_ps_dan01_all+1), cluster_rows=F, cluster_cols=F,scale = "row",fontsize=2)

## ODCs:
sn_atlas_processed_tmp<-subset(sn_atlas_processed,subset=CellType=="ODC")
sn_atlas_processed_tmp <- NormalizeData(sn_atlas_processed_tmp)
sn_atlas_processed_tmp <- FindVariableFeatures(sn_atlas_processed_tmp, selection.method = "vst", nfeatures = 2000)
sn_atlas_processed_tmp <- ScaleData(sn_atlas_processed_tmp)
sn_atlas_processed_tmp <- RunPCA(sn_atlas_processed_tmp, features = VariableFeatures(object = sn_atlas_processed_tmp))
sn_atlas_processed_tmp <- FindNeighbors(sn_atlas_processed_tmp, dims = 1:30)
sn_atlas_processed_tmp <- FindClusters(sn_atlas_processed_tmp, resolution = 0.6)
sn_atlas_processed_tmp <- RunUMAP(sn_atlas_processed_tmp, dims = 1:30)
sce_odc <- as.SingleCellExperiment(sn_atlas_processed_tmp, assay = "RNA")
sce_odc <- slingshot(sce_odc, reducedDim = 'UMAP', clusterLabels = 'Disease')


meta_toplot_odc<-data.frame(ps=sce_odc$slingPseudotime_1,disease=sce_odc$Disease,cellt=sce_odc$CellSubType)
meta_toplot_odc<-meta_toplot_odc[meta_toplot_odc$cellt=="ODC_2",]
ggplot(meta_toplot_odc, aes(x=ps, fill=disease)) +
    geom_density(alpha=.5) + theme_bw()


# characterisation of the genes changing expression along pseudotime

genes_odc2<-as.matrix(GetAssayData(sn_atlas_processed_tmp, slot = "counts"))
sce_odc2_deg <- fitGAM(genes_odc2,sds=SlingshotDataSet(sce_odc2))
odc2_deg_tradeseq <- associationTest(sce_odc2_deg)
odc2_deg_tradeseqStartEnd <- startVsEndTest(sce_odc2_deg)

save(sce_odc,sce_odc2_deg,odc2_deg_tradeseqStartEnd,file="results/pseudotime/ps_odc.RData")



