library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(igraph)
library(here)

# Load DEGs
load(here("results/DEG_tables/DEG_results.RData"))

# Function to perform GO enrichment for a list of genes
go_enrichment <- function(gene_list, ont = "BP") {
  gene_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  enrichGO(gene = gene_ids$ENTREZID,
           OrgDb = org.Hs.eg.db,
           ont = ont,
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.1)
}

# Loop over each cluster's DEGs
go_results_degs <- list()
for (clust in names(deg_results)) {
  genes <- rownames(deg_results[[clust]])[deg_results[[clust]]$p_val_adj < 0.05]
  if(length(genes) > 0) {
    go_results_degs[[clust]] <- go_enrichment(genes)
  }
}



# Cell type-specific PPI networks
# to be done for DaN_3, ODC_2 and OPC (where PD genetic risk converged)

ppi_net<-read.delim("ppi_net_new2023short",h=F)
ppi_net_dan3top1000<-ppi_net[ppi_net[,1]%in%rownames(df_w_test_DaN[oder(-df_w_test_DaN$DaN_3),])[1:1000] & ppi_net[,2]%in%rownames(df_w_test_DaN[oder(-df_w_test_DaN$DaN_3),])[1:1000] ,]
ppi_net_dan3top1000_genes<-unique(c(as.character(ppi_net_dan3top1000[,1]),as.character(ppi_net_dan3top1000[,2])))

net_1 <- graph_from_data_frame(d= ppi_net_dan3top1000, vertices= ppi_net_dan3top1000_genes, directed=F)
cl_net_1<-cluster_louvain(net_1)
memb_net_1<-as.matrix(membership(cl_net_1))
memb_net_dan3top1000ppi_df<-data.frame(id=rownames(memb_net_1),module=memb_net_1)
ppi_results <- memb_net_dan3top1000ppi_df

# Loop over each gene module
go_results_ppi_network <- list()
for (module in unique(ppi_results$module)) {
  genes <- ppi_results[ppi_results$module==module,1]
  if(length(genes) > 0) {
    go_results_ppi_network[[module]] <- go_enrichment(genes)
  }
}

save(go_results_degs, go_results_ppi_network, file = here("results/GO_PPI/GO_PPI_results.RData"))
