library(Seurat)
library(dplyr)
library(here)

# Load subclustered object
load(here("data/processed/DA_subclusters.RData"))

### Find cell type specific gene markers 

```{r gene markers, echo=TRUE}

# Method 1 (used in the manuscript):
# to be done for each cell type and subtype in controls

# level 1:

list_w_stat<-c()
df_w_test<-c()
signed_rank = function(x) sign(x) * rank(abs(x))
metadataALL<-as.data.frame(sn_atlas_processed@meta.data)
metadataALL_CTR<-metadataALL[metadataALL$Disease=="CTR",]
data<-as.matrix(GetAssayData(sn_atlas_processed, slot = "counts"))
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
    return(val$coefficients[2,3]) ## w-statisc value
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

# level 2:

list_w_stat<-c()
df_w_test<-c()
signed_rank = function(x) sign(x) * rank(abs(x))
metadataALL<-as.data.frame(sn_atlas_processed@meta.data)
metadataALL_CTR<-metadataALL[metadataALL$Disease=="CTR",]
metadataALL_CTR<-metadataALL_CTR[metadataALL_CTR$CellType=="ODC",]
data<-as.matrix(GetAssayData(sn_atlas_processed, slot = "counts"))
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
    return(val$coefficients[2,3]) ## w-statisc value
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
write.table(df_w_test,"Wtest_ODCs",quote=F,sep="\t")





# Differential expression
deg_results <- list()
for (clust in levels(Idents(DA_cells))) {
  deg_results[[clust]] <- FindMarkers(DA_cells, ident.1 = clust, min.pct = 0.25)
}

# Save DEG results
save(deg_results, file = here("results/DEG_tables/DA_DEG_results.RData"))
