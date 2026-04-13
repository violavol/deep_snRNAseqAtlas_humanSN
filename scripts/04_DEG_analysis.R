library(Seurat)
library(dplyr)
library(here)

# Load subclustered object
load(here("data/processed/sn_atlas_annotated_subtype.RData"))

# Find cell type specific gene markers for MAGMA analysis

```{r gene markers, echo=TRUE}


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
save(df_w_test_level1, df_w_test_level2_DaN, file = here("results/DEG_tables/DEG_results.RData"))
