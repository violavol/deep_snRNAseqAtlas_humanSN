---

## File 2 — `R/cellchat_utils.R`

```r

prepare_cellchat_object <- function(
  seurat_obj,
  disease_value,
  disease_col,
  celltype_col,
  downsample_n,
  db = CellChatDB.human
) {
  obj_sub <- subset(seurat_obj, subset = .data[[disease_col]] == disease_value)
  Idents(obj_sub) <- obj_sub[[celltype_col]][, 1]
  obj_sub <- subset(obj_sub, downsample = downsample_n)

  expr <- GetAssayData(obj_sub, assay = "RNA", slot = "data")
  labels <- obj_sub[[celltype_col]][, 1]
  names(labels) <- colnames(expr)

  meta <- data.frame(group = labels, row.names = names(labels))
  cellchat_obj <- createCellChat(object = expr, meta = meta, group.by = "group")
  cellchat_obj@DB <- db

  cellchat_obj
}

run_cellchat_pipeline <- function(cellchat_obj, min_cells = 10, raw_use = TRUE, population_size = FALSE) {
  cellchat_obj <- subsetData(cellchat_obj)
  cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
  cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
  cellchat_obj <- computeCommunProb(cellchat_obj, raw.use = raw_use, population.size = population_size)
  cellchat_obj <- filterCommunication(cellchat_obj, min.cells = min_cells)
  cellchat_obj <- computeCommunProbPathway(cellchat_obj)
  cellchat_obj <- aggregateNet(cellchat_obj)
  cellchat_obj <- netAnalysis_computeCentrality(cellchat_obj, slot.name = "netP")

  cellchat_obj
}

export_cellchat_tables <- function(cellchat_obj, prefix, output_dir) {
  df_net <- subsetCommunication(cellchat_obj)
  df_netP <- subsetCommunication(cellchat_obj, slot.name = "netP")

  write.table(
    df_net,
    file = file.path(output_dir, paste0(prefix, "_communication_network.tsv")),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )

  write.table(
    df_netP,
    file = file.path(output_dir, paste0(prefix, "_signaling_pathways.tsv")),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )

  list(df_net = df_net, df_netP = df_netP)
}

plot_signaling_role_heatmaps <- function(cellchat_obj) {
  ht_out <- netAnalysis_signalingRole_heatmap(
    cellchat_obj,
    pattern = "outgoing",
    font.size = 4,
    height = 20,
    width = 6
  )

  ht_in <- netAnalysis_signalingRole_heatmap(
    cellchat_obj,
    pattern = "incoming",
    font.size = 4,
    height = 20,
    width = 6
  )

  list(outgoing = ht_out, incoming = ht_in)
}

merge_cellchat_objects <- function(cellchat_ctrl, cellchat_pd) {
  common_groups <- levels(cellchat_ctrl@idents)
  cellchat_pd <- liftCellChat(cellchat_pd, common_groups)

  mergeCellChat(
    list(CTR = cellchat_ctrl, PD = cellchat_pd),
    add.names = c("CTR", "PD")
  )
}

run_interaction_enrichment_test <- function(
  cellcomm,
  dan_genes,
  pd_risk_genes,
  source_cell,
  target_cell,
  n_perm = 10000
) {
  interact_1 <- cellcomm[
    cellcomm$target == target_cell &
      cellcomm$source == source_cell &
      cellcomm$receptor %in% dan_genes,
  ]

  interact_2 <- cellcomm[
    cellcomm$source == target_cell &
      cellcomm$target == source_cell &
      cellcomm$ligand %in% dan_genes,
  ]

  interact_bg_1 <- cellcomm[
    cellcomm$target == target_cell &
      cellcomm$source == source_cell,
  ]

  interact_bg_2 <- cellcomm[
    cellcomm$source == target_cell &
      cellcomm$target == source_cell,
  ]

  observed_1 <- interact_1[interact_1$ligand %in% pd_risk_genes, ]
  observed_2 <- interact_2[interact_2$receptor %in% pd_risk_genes, ]
  observed_hits <- unique(c(observed_1$ligand, observed_2$receptor))

  if (length(observed_hits) == 0) {
    return(list(
      p_value = 1,
      n_observed = 0,
      observed_hits = character(0)
    ))
  }

  n_more_extreme <- 0

  for (i in seq_len(n_perm)) {
    sampled_risk <- sample(pd_risk_genes, length(pd_risk_genes))
    sampled_ligands <- sample(interact_bg_1$ligand, length(unique(interact_1$ligand)))
    sampled_receptors <- sample(interact_bg_2$receptor, length(unique(interact_2$receptor)))

    perm_hits_1 <- sampled_ligands[sampled_ligands %in% sampled_risk]
    perm_hits_2 <- sampled_receptors[sampled_receptors %in% sampled_risk]
    perm_hits <- unique(c(perm_hits_1, perm_hits_2))

    if (length(perm_hits) > length(observed_hits)) {
      n_more_extreme <- n_more_extreme + 1
    }
  }

  list(
    p_value = n_more_extreme / n_perm,
    n_observed = length(observed_hits),
    observed_hits = observed_hits
  )
}
