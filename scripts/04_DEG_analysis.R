library(Seurat)
library(dplyr)
library(here)

# Load subclustered object
load(here("data/processed/sn_atlas_annotated_subtype.RData"))

# Find DEGs by subtype
sn_combined$clust.disease <- paste(sn_combined$CellSubType, sn_combined$Disease, sep = "_")
Idents(sn_combined) <- "clust.disease"

sn_combined <- subset(sn_combined,subset=Disease=="CTR" | Disease=="PD_B5-6")
deg_results <- list()
for (clust in levels(sn_combined)) {
    markers <- FindMarkers(
    sn_combined,
    ident.1 = clust,
    group.by = "Disease"   # CTR vs PD_B5_6
  )
  deg_results[[clust]] <- markers
}

# Save  results
save(deg_results, df_w_test_level2_DaN, file = here("results/DEG_tables/DEG_results.RData"))



# Find cell type specific gene markers for MAGMA analysis
# to be done for each cell type and subtype in control samples

# level 1:

list_w_stat<-c()
df_w_test<-c()
signed_rank = function(x) sign(x) * rank(abs(x))
metadataALL<-as.data.frame(sn_combined@meta.data)
metadataALL_CTR<-metadataALL[metadataALL$Disease=="CTR",]
data<-as.matrix(GetAssayData(sn_combined, slot = "counts"))
data<-data[,colnames(data)%in%rownames(metadataALL_CTR)]
list_cell <- as.array(unique(metadataALL_CTR$CellType))
list_gene <- row.names(data)

compute_w_val<-function(gene_ref,cell_ref,matrix_expr,info_cell) {
    n_cell<- ncol(data)
    X<- matrix(-1, nrow = n_cell, ncol = 1)
    X[which(info_cell$CellSubType==cell_ref),1]<-1
    Y=t(data[gene_ref,])
    df<-data.frame(y=t(Y),x=X)
    colnames(df)<-c("y","x")
    linearMod <- lm(signed_rank(y) ~ x, data=df)
    val<-summary(linearMod)
    list_val<-val$coefficients[2,3]
    return(val$coefficients[2,3]) ## statistic
}

for (i in c(1:length(list_cell))) {
    cell_ref=list_cell[i]
    print("W test")
    print(list_cell[i])
    list_w_stat<-sapply(list_gene,compute_w_val,cell_ref=list_cell[i],matrix_expr=data,info_cell=metadataALL_CTR)
    if(!exists("df_w_test")) {
        df_w_test<-data.frame(list_w_stat)
    } else {
        df_w_test<-cbind(df_w_test,list_w_stat)
    }
}
colnames(df_w_test) <- list_cell
write.table(df_w_test,"Wtest_level1",quote=F,sep="\t")
df_w_test_level1 <- df_w_test

# level 2 (DaN subtype):

list_w_stat<-c()
df_w_test<-c()
signed_rank = function(x) sign(x) * rank(abs(x))
metadataALL<-as.data.frame(sn_combined@meta.data)
metadataALL_CTR<-metadataALL[metadataALL$Disease=="CTR",]
metadataALL_CTR<-metadataALL_CTR[metadataALL_CTR$CellType=="DaN",]
data<-as.matrix(GetAssayData(sn_combined, slot = "counts"))
data<-data[,colnames(data)%in%rownames(metadataALL_CTR)]
list_cell <- as.array(unique(metadataALL_CTR$CellSubType))
list_gene <- row.names(data)

compute_w_val<-function(gene_ref,cell_ref,matrix_expr,info_cell) {
    n_cell<- ncol(data)
    X<- matrix(-1, nrow = n_cell, ncol = 1)
    X[which(info_cell$CellSubType==cell_ref),1]<-1
    Y=t(data[gene_ref,])
    df<-data.frame(y=t(Y),x=X)
    colnames(df)<-c("y","x")
    linearMod <- lm(signed_rank(y) ~ x, data=df)
    val<-summary(linearMod)
    list_val<-val$coefficients[2,3]
    return(val$coefficients[2,3]) ## statistic
}

for (i in c(1:length(list_cell))) {
    cell_ref=list_cell[i]
    print("W test")
    print(list_cell[i])
    list_w_stat<-sapply(list_gene,compute_w_val,cell_ref=list_cell[i],matrix_expr=data,info_cell=metadataALL_CTR)
    if(!exists("df_w_test")) {
        df_w_test<-data.frame(list_w_stat)
    } else {
        df_w_test<-cbind(df_w_test,list_w_stat)
    }
}
colnames(df_w_test) <- list_cell
write.table(df_w_test,"Wtest_DaN_subtypes",quote=F,sep="\t")
df_w_test_level2_DaN <- df_w_test


# Save  results
save(df_w_test_level1, df_w_test_level2_DaN, file = here("results/DEG_tables/DEG_results_MAGMA.RData"))



# Prepare for MAGMA analysis 
# level 1:

df_w_test_2p<-data.frame(cellt=c(rep("ODC",500),rep("Astrocyte",500),rep("DaN",500),rep("GABA",500),rep("OPC",500),rep("Microglia",500),rep("Tcell",500)),id=c(rownames(df_w_test[order(-df_w_test$ODC),])[1:500],rownames(df_w_test[order(-df_w_test$Astrocyte),])[1:500],rownames(df_w_test[order(-df_w_test$DaN),])[1:500],rownames(df_w_test[order(-df_w_test$GABA),])[1:500],rownames(df_w_test[order(-df_w_test$OPC),])[1:500],rownames(df_w_test[order(-df_w_test$Microglia),])[1:500],rownames(df_w_test[order(-df_w_test$Tcell),])[1:500]))

NCBI37.3.gene.loc<-read.delim("NCBI37.3.gene.loc",h=F)
names(NCBI37.3.gene.loc)[6]<-"id"
df<-merge(df_w_test_2p[,1:2],NCBI37.3.gene.loc,by=c("id"))
list_gene<-c()
for(m in unique(df$cellt)){
 val<-c(m,as.character(df[df$cellt==m,3]))
 val<-paste(val,collapse =" ")
 list_gene<-c(list_gene,val)
}
write.table(data.frame(list_gene),file="input_magma_level1_2p",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
