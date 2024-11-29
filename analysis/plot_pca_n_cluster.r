# Load library ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(here)
  library(sessioninfo)
})

# Load data ----
spe_joint <- readRDS(
  here("processed_data", "spe_joint.rds")
)


# Make plots ----
## Visium (2000 HVG) with Joint PCA ----

# TODO: delete the commented out session
# spe_joint$cluster_HVG_label <- c(
#   spe_visium$cluster_HVG_label,
#   rep(NA, ncol(spe_merf))
# ) |> factor()

spe_joint[, spe_joint$sample_id == "visium"] |>
  plotReducedDim(
    ncomponents = 2,
    dimred = "PCA",
    color_by = "cluster_HVG_label" # ,
    # point_size = 0.5
  ) +
  scale_x_continuous(
    limits = c(-10, 10)
  )


## Visium (overlap genes) with Joint PCA ----

# spe_joint$cluster_400_label <- c(
#   spe_visium$cluster_400_label,
#   rep(NA, ncol(spe_merf))
# ) |> factor()

spe_joint[, spe_joint$sample_id == "visium"] |>
  plotReducedDim(
    ncomponents = 2,
    dimred = "PCA",
    color_by = "cluster_400_label",
    point_size = 0.5
  )


## Merfish (overlap genes) with Joint PCA ----
# spe_joint$merf_label <- c(
#   rep(NA, ncol(spe_visium)),
#   spe_merf$merf_label
# ) |> factor()

p_merf <- spe_joint[, spe_joint$sample_id == "merfish"] |>
  plotReducedDim(
    ncomponents = 2,
    dimred = "PCA",
    color_by = "merf_label",
    point_size = 0.5
  )


# Adjust color ----

## Find aligned color between HVG cluster and 2000 HVG cluster ----
# TODO: Heatmap:
# https://bioconductor.org/books/3.20/OSCA.basic/clustering.html#adjusting-the-parameters



# Create Planel ----


# Save plot output ----

# Session Info ----
session_diff()
