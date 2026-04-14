### Coexpression network analysis

```{r coexpr, echo=TRUE}
library(Seurat)
library(BigScale2)


load(here("data/processed/sn_atlas_annotated_subtype.RData"))
sn_atlas_DAN3 <- subset(sn_combined,subset=CellSubType=="DaN_3")

data_dan3<-as.matrix(GetAssayData(sn_atlas_DAN3, slot = "counts"))
results.dan3=compute.network(expr.data = data_dan3,gene.names = rownames(data_dan3),clustering = "direct")

m<-data.frame(results.dan3$correlations)
coexp_dan3<-data.frame(row=rownames(m)[row(m)[upper.tri(m)]], 
                       col=colnames(m)[col(m)[upper.tri(m)]], 
                       corr=m[upper.tri(m)])


coexp_dan3_08<-coexp_dan3[abs(coexp_dan3$corr)>0.8,]

coexp_dan3_08_genes<-unique(c(as.character(coexp_dan3_08[,1]),as.character(coexp_dan3_08[,2])))
net_dan3 <- graph_from_data_frame(d= coexp_dan3_08, vertices= coexp_dan3_08_genes, directed=F)
cl_net<-cluster_louvain(net_dan3)
memb_net_dan3<-as.matrix(membership(cl_net))
memb_net_dan3_df<-data.frame(id=rownames(memb_net_dan3),module=memb_net_dan3)
