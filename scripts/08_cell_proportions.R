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
df_prop <- meta %>%
  count(CellSubType, Disease) %>%          # count cells
  group_by(CellSubType) %>%
  mutate(prop = n / sum(n)) %>%       # compute proportions
  ungroup()


prop_plot <- ggplot(df_prop, aes(x = celltype, y = prop, fill = group)) +
  geom_col(width = 0.9) +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_manual(
    values = c("CTR" = "mediumorchid", "PD_B5_6" = "skyblue2")
  ) +
  labs(x = NULL, y = "prop", fill = NULL) +
  theme_gray() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )



# Save figure
ggsave(here("figures/Figure_3/CellType_Proportions.pdf"), prop_plot, width = 6, height = 4)
save(df_prop, file = here("results/CellType_Proportions.RData"))
