library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(igraph)
library(here)

# Load DEGs
load(here("results/DEG_tables/DA_DEG_results.RData"))

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
go_results <- list()
for (clust in names(deg_results)) {
  genes <- rownames(deg_results[[clust]])[deg_results[[clust]]$p_val_adj < 0.05]
  if(length(genes) > 0) {
    go_results[[clust]] <- go_enrichment(genes)
  }
}

# Optional: construct PPI network for top DEGs using STRINGdb (or igraph)
# Here we create a simple example with igraph
ppi_edges <- data.frame(from=c(), to=c())
for (clust in names(go_results)) {
  top_genes <- head(rownames(deg_results[[clust]]), 20)
  edges <- expand.grid(top_genes, top_genes)
  edges <- edges[edges$Var1 != edges$Var2, ]
  ppi_edges <- rbind(ppi_edges, edges)
}

ppi_network <- graph_from_data_frame(ppi_edges, directed = FALSE)
save(go_results, ppi_network, file = here("results/GO_PPI/GO_PPI_results.RData"))
