# Load library ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(here)
  library(scran)
  library(scater)
  library(sessioninfo)
  library(pheatmap)
  library(tidyverse)
})

# Load data ----
spe_joint <- readRDS(
  here("processed_data", "spe_joint.rds")
)


set.seed(1100101001)
spe_joint <- runUMAP(spe_joint, dimred = "PCA")
# plotReducedDim(spe_joint, dimred="UMAP", colour_by="level1class")



# Make plots ----
## Color match ----
matched_df <- tribble(
  ~overlap, ~HVG, ~overlap_color, ~HVG_color,
  2, 1, "#DE659F", "#EBADC9",
  2, 3, "#DE659F", "#DE659F",
  9, 5, "#B64F31", "#B64F31",
  6, 4, "#603B8F", "#603B8F",
  3, 12, "#F8D949", "#F8D949",
  5, 16, "#D75B42", "#D75B42",
  5, 17, "#D75B42", "#E18D7E",
  7, 9, "#181818", "#181818",
  8, 15, "#5FAE77", "#5FAE77",
  1, 11, "#D17739", "#EED1AE",
  1, 13, "#D17739", "#E0A261",
  1, 7, "#D17739", "#D17739",
  4, 10, "#4A7DB4", "#173754",
  4, 2, "#4A7DB4", "#255989",
  4, 6, "#4A7DB4", "#306FAB",
  4, 8, "#4A7DB4", "#508EBC",
) |> mutate(
  overlap_lab = paste("Overlap", overlap),
  HVG_lab = paste("HVG", HVG),
)



## Heatmap for label annotation & Color selection ----
# Using a large pseudo-count for a smoother color transition
spe_visium <- spe_joint[, spe_joint$sample_id == "visium"]

# between 0 and 1 cell in each 'tab'.
tab <- table(
  paste("HVG", spe_visium$cluster_HVG_label),
  paste("Overlap", spe_visium$cluster_400_label)
)

anno_colors <- list(
  overlap_lab = matched_df |>
    group_by(overlap_lab) |> slice_head(n = 1) |> ungroup() |>
    with(setNames(overlap_color, nm = overlap_lab)),
  HVG_lab = with(matched_df, setNames(HVG_color, nm = HVG_lab))
)

anno_row_df <- matched_df |>
  group_by(overlap_lab) |>
  slice_head(n = 1) |>
  ungroup() |>
  pull(overlap_lab) |>
  setNames(nm = _) |>
  data.frame()
colnames(anno_row_df) <- "overlap_lab"

anno_col_df <- matched_df |>
  pull(HVG_lab) |>
  setNames(nm = _) |>
  data.frame()

colnames(anno_col_df) <- "HVG_lab"

pdf(
  here(
    "plots/case_visium_n_merf",
    "heatmap_match_cluster_labels.pdf"
  ),
  height = 1.6 * 2,
  width = 2.8 * 2
)
pheatmap(
  log10(tab + 10) |> t(), # Row: overlap, #column HVG
  # main = "Infomap vs Walktrap",
  color = viridis::viridis(100),
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  annotation_row = anno_row_df,
  annotation_col = anno_col_df,
  annotation_colors = anno_colors,
  annotation_legend = FALSE,
  annotation_names_row = FALSE, annotation_names_col = FALSE,
  annotation_row_side = "bottom", annotation_col_side = "right",
  treeheight_col = 0 # Hide the column dendrogram
)
dev.off()


## Visium (2000 HVG) with Joint PCA ----
# TODO: delete the commented out session
# spe_joint$cluster_HVG_label <- c(
#   spe_visium$cluster_HVG_label,
#   rep(NA, ncol(spe_merf))
# ) |> factor()

umap_hvg <- spe_joint[, spe_joint$sample_id == "visium"] |>
  plotReducedDim(
    ncomponents = 2,
    dimred = "UMAP",
    color_by = "cluster_HVG_label",
    point_size = 0.5
  ) +
  # scale_x_continuous(
  #   limits = c(-10, 10)
  # ) +
  scale_color_manual(
    values = setNames(
      matched_df$HVG_color,
      matched_df$HVG |> as.character()
    ),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(-1, NA), ylim = c(7.5, NA)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )


## Visium (overlap genes) with Joint PCA ----

# spe_joint$cluster_400_label <- c(
#   spe_visium$cluster_400_label,
#   rep(NA, ncol(spe_merf))
# ) |> factor()

umap_overlap <- spe_joint[, spe_joint$sample_id == "visium"] |>
  plotReducedDim(
    ncomponents = 2,
    dimred = "UMAP",
    color_by = "cluster_400_label",
    point_size = 0.5
  ) +
  # scale_x_continuous(
  #   limits = c(-10, 10)
  # ) +
  scale_color_manual(
    values = setNames(
      matched_df$overlap_color,
      matched_df$overlap |> as.character()
    ),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(-1, NA), ylim = c(7.5, NA)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

pdf(
  here(
    "plots/case_visium_n_merf",
    "umap_visium_cluster_diff.pdf"
  ),
  heigh = 2,
  width = 4.5
)
ggpubr::ggarrange(
  umap_overlap,
  umap_hvg
)
dev.off()

## Merfish (overlap genes) with Joint PCA ----
# spe_joint$merf_label <- c(
#   rep(NA, ncol(spe_visium)),
#   spe_merf$merf_label
# ) |> factor()

p_merf <- spe_joint[, spe_joint$sample_id == "merfish"] |>
  plotReducedDim(
    ncomponents = 2,
    dimred = "UMAP",
    color_by = "merf_label",
    point_size = 0.5
  )

## Colored by tech ----
spe_joint |>
  plotReducedDim(
    ncomponents = 2,
    dimred = "UMAP",
    color_by = "sample_id",
    point_size = 0.2
  ) +
  scale_color_manual(
    values = c(
      "merfish" = "#5e3c99", "visium" = "#fdb863"
    ), name = "Data",
    labels = c("MERFISH", "Visium"),
    guide = "none"
  ) #+
# theme(legend.text = element_text(size = 11)) +
# guides(color = guide_legend(override.aes = list(size = 3))) # Adjust the size as needed

# Save the plot to a PDF file
ggsave(
  filename = here(
    "plots/case_visium_n_merf",
    "output_figure.pdf"
  ),
  plot = last_plot(),
  device = "pdf",
  width = 3,
  height = 2,
  units = "in"
)
# theme_bw(base_size = 12)
# # scale_x_continuous(
# #   limits = c(-10, 10)
# # ) +
# scale_color_manual(
#   values = setNames(
#     matched_df$overlap_color,
#     matched_df$overlap |> as.character()
#   )
# )


# Spatial Plots ----
library(escheR)

pdf(
  here("plots/case_visium_n_merf",
  "spatial_plot_visium_cluster_overlap.pdf"),
  width = 4.2,
  height = 4.2
)
make_escheR(spe_visium) |>
  add_fill(var = "cluster_400_label", point_size = 2.2) +
  scale_fill_manual(
    values = setNames(
      matched_df$overlap_color,
      matched_df$overlap |> as.character()
    ),
    guide = "none"
  )
dev.off()

pdf(
  here("plots/case_visium_n_merf",
  "spatial_plot_visium_cluster_HVG.pdf"),
  width = 4.2,
  height = 4.2
)
make_escheR(spe_visium) |>
  add_fill(var = "cluster_HVG_label") +
  scale_fill_manual(
    values = setNames(
      matched_df$HVG_color,
      matched_df$HVG |> as.character()
    ),
    guide = "none"
  )
dev.off()
# Adjust color ----

## Find aligned color between HVG cluster and 2000 HVG cluster ----
# TODO: Heatmap:
# https://bioconductor.org/books/3.20/OSCA.basic/clustering.html#adjusting-the-parameters



# Create Planel ----


# Save plot output ----

# Session Info ----
session_diff()
