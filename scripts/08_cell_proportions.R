# 08_cell_proportions.R
# Differential cell-type proportions between CTR and PD_B5-6
# using propeller and scProportionTest.
# Project: deep_snRNAseqAtlas_humanSN — Volpato et al. (2026)

library(Seurat)
library(dplyr)
library(ggplot2)
library(here)
library(speckle)
library(limma)
library(scProportionTest)
source(here("scripts/utils.R"))

# ── Parameters ────────────────────────────────────────────────────────────────
CTRL_LABEL <- "CTR"
PD_LABEL   <- "PD_B5-6"
COLOR_MAP  <- c("CTR" = "mediumorchid", "PD_B5-6" = "skyblue2")

# ── Load data ─────────────────────────────────────────────────────────────────
load(here("data/processed/sn_atlas_annotated_subtype.RData"))

# ── propeller ─────────────────────────────────────────────────────────────────
propeller_res <- propeller(
  clusters = sn_combined$CellSubType,
  sample   = sn_combined$Sample_v2,
  group    = sn_combined$Disease
)

# ── scProportionTest ──────────────────────────────────────────────────────────
prop_test <- sc_utils(sn_combined)
prop_test <- permutation_test(
  prop_test,
  cluster_identity  = "CellSubType",
  sample_1          = CTRL_LABEL,
  sample_2          = PD_LABEL,
  sample_identity   = "Disease"
)

# ── Proportion data-frame ─────────────────────────────────────────────────────
meta <- sn_combined@meta.data

df_prop <- meta %>%
  count(CellSubType, Disease) %>%
  group_by(CellSubType) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# ── Plot ──────────────────────────────────────────────────────────────────────
# Fix variable name bug from original: use CellSubType/Disease (not celltype/group)
prop_plot <- ggplot(df_prop, aes(x = CellSubType, y = prop, fill = Disease)) +
  geom_col(width = 0.9) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = COLOR_MAP) +
  labs(x = NULL, y = "Proportion (%)", fill = "Disease") +
  theme_gray() +
  theme(
    axis.text.x    = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# ── Save ──────────────────────────────────────────────────────────────────────
ensure_dirs(here("figures/Figure_3"), here("results"))
save_plot(prop_plot, here("figures/Figure_3/CellType_Proportions.pdf"),
          width = 6, height = 4)

save(df_prop, propeller_res, prop_test,
     file = here("results/CellType_Proportions.RData"))
message("Cell proportion results saved.")

log_session("08_cell_proportions")
