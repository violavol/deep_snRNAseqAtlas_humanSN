library(Seurat)
library(dplyr)
library(ggplot2)
library(here)

# Load integrated annotated object
load(here("data/processed/sn_atlas_annotated.RData"))

# Compute cell type proportions per disease group
cell_counts <- table(sn_combined$Disease, Idents(sn_combined))
cell_props <- prop.table(cell_counts, margin = 1)
cell_props_df <- as.data.frame(cell_props)
colnames(cell_props_df) <- c("Disease", "CellType", "Proportion")

# Plot cell type proportions
prop_plot <- ggplot(cell_props_df, aes(x=Disease, y=Proportion, fill=CellType)) +
  geom_bar(stat="identity", position="stack") +
  theme_minimal() +
  labs(title="Cell Type Proportions per Disease Group")

# Save figure
ggsave(here("figures/Figure_3/CellType_Proportions.pdf"), prop_plot, width = 6, height = 4)
save(cell_props_df, file = here("results/CellType_Proportions.RData"))
