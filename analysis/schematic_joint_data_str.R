# Load library ----
library(here)
library(SpatialExperiment)
library(pheatmap)
library(sessioninfo)

# Load Data ----
spe_visium <- readRDS(
  here("processed_data", "spe_visium.rds")
)

spe_merf <- readRDS(
  here("processed_data", "spe_merf.rds")
)

# Extract the count matrices
counts_visium <- assay(spe_visium, "counts") |> data.matrix()
counts_merf <- assay(spe_merf, "counts") |> data.matrix()


merf_unique_gene <- setdiff(rownames(spe_merf), rownames(spe_visium))
overlap_gene <- intersect(rownames(spe_merf), rownames(spe_visium))
visium_unique_gene <- setdiff(rownames(spe_visium), rownames(spe_merf))

set.seed(20241129)
merf_unique_gene_sampled <- sample(merf_unique_gene, size = 0.01 * length(merf_unique_gene))
overlap_gene_sampled <- sample(overlap_gene, size = 0.01 * length(overlap_gene))
visium_unique_gene_sampled <- sample(visium_unique_gene, size = 0.01 * length(visium_unique_gene))

merf_height <- ceiling(0.01 * ncol(spe_merf))
visium_height <- ceiling(0.01 * ncol(spe_visium))

merf_df <- cbind(
  matrix(-1,
    ncol = length(merf_unique_gene_sampled) + length(overlap_gene_sampled),
    nrow = merf_height
  ) +
    rnorm(
      n = merf_height * (length(merf_unique_gene_sampled) + length(overlap_gene_sampled)),
      sd = 0.1
    ),
  matrix(NA,
    ncol = length(visium_unique_gene_sampled),
    nrow = merf_height
  )
)

visium_df <- cbind(
  matrix(NA,
    ncol = length(merf_unique_gene_sampled),
    nrow = visium_height
  ),
  matrix(0.5,
    ncol = length(overlap_gene_sampled),
    nrow = visium_height
  ) +
    rnorm(
      n = visium_height * length(overlap_gene_sampled),
      sd = 0.005
    ),
  matrix(1,
    ncol = length(visium_unique_gene_sampled),
    nrow = visium_height
  ) +
    rnorm(
      n = visium_height * length(visium_unique_gene_sampled),
      sd = 0.1
    )
)

joint_df <- rbind(
  merf_df,
  visium_df
)

# Make data heatmap ----
pheatmap(
  joint_df,
  color = colorRampPalette(c("#5e3c99", "white", "#e66101"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  na_col = "#ECECEC",
  legend = FALSE
) |> print()

# Session Info ----
session_info()
