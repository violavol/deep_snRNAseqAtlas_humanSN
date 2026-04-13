library(Seurat)
library(dplyr)
library(ggplot2)
library(here)
library(speckle)
library(limma)
library(ggplot2)
library("scProportionTest")

# Load integrated annotated object
load(here("data/processed/sn_atlas_annotated_subtype.RData"))

# Compute cell type proportions using propeller
propeller(clusters = sn_combined$CellSubType, sample = sn_combined$Sample_v2, 
group = sn_combined$Disease)

# Compute cell type proportions using scProportionTest
prop_test <- sc_utils(sn_combined)
prop_test <- permutation_test(
	prop_test, cluster_identity = "CellSubType",
	sample_1 = "CTR", sample_2 = "PD_B5-6",
	sample_identity = "Disease"
)


# Plot cell type proportions
meta <- sn_combined@meta.data

df_counts <- meta %>%
    group_by(donor_id, disease, seurat_clusters) %>%
    summarise(n = n(), .groups = "drop")
df_prop <- df_counts %>%
    group_by(donor_id) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()

ggplot(df_prop, aes(x = disease, y = prop, fill = disease)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~ seurat_clusters, scales = "free_y") +
    theme_classic() +
    ylab("Cell type proportion") +
    xlab("Disease")




prop_plot <- ggplot(cell_props_df, aes(x=Disease, y=Proportion, fill=CellType)) +
  geom_bar(stat="identity", position="stack") +
  theme_minimal() +
  labs(title="Cell Type Proportions per Disease Group")

# Save figure
ggsave(here("figures/Figure_3/CellType_Proportions.pdf"), prop_plot, width = 6, height = 4)
save(cell_props_df, file = here("results/CellType_Proportions.RData"))
