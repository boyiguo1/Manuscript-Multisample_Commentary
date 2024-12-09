# Load library ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(here)
  library(scran)
  library(scater)
  library(sessioninfo)
  library(pheatmap)
  library(tidyverse)
  library(escheR)
})

# Create SPE object ----

# Merfish data ----
# NOTE: we used the results from Clifton et al. (2023)
merf_csv <- readr::read_csv(here("raw_data/STalign_S2R3_to_Visium.csv"))
colnames(merf_csv)[1] <- "key"
count_assay <- merf_csv[, 7:ncol(merf_csv)] |> t() # Gene-by-cell
cd <- merf_csv[, 1:6]

raw_merf <-
  pdf(
    here(
      "plots/case_visium_n_merf",
      "spatial_plot_raw_merfish.pdf"
    ),
    height = 1 * 3, width = 1.45 * 3
  )
SpatialExperiment(
  assay = list("counts" = count_assay),
  sample_id = "merfish",
  colData = cd,
  spatialCoordsNames = c("center_x", "center_y")
) |>
  make_escheR() |>
  add_fill("sample_id", point_size = 0.5) +
  scale_fill_manual(
    values = c("merfish" = "#5e3c99"), guide = "none"
  ) + scale_y_continuous()
dev.off()


pdf(
  here(
    "plots/case_visium_n_merf",
    "spatial_plot_STaligned_merfish.pdf"
  ),
  height = 1 * 3, width = 1.45 * 3
)
SpatialExperiment(
  assay = list("counts" = count_assay),
  sample_id = "merfish",
  colData = cd,
  spatialCoordsNames = c("STalign_x", "STalign_y")
) |>
  make_escheR() |>
  add_fill("sample_id", point_size = 0.5) +
  scale_fill_manual(
    values = c("merfish" = "#5e3c99"), guide = "none"
  )
dev.off()


# rowData(raw_merf)$symbol <- rownames(raw_merf)
# # Create mouse ensembl id in row data
# ensmbl_vec <- rowData(spe_visium)$ensmbl[match(rowData(raw_merf)$symbol, rowData(spe_visium)$symbol)]
# rowData(raw_merf)$ensmbl <- ensmbl_vec
# # Non-matching ensmbl use gene symbol
# rowData(raw_merf)$ensmbl[which(is.na(ensmbl_vec))] <- rowData(raw_merf)$symbol[which(is.na(ensmbl_vec))]

# # Use ensmbl id as row name
# rownames(raw_merf) <- rowData(raw_merf)$ensmbl

# # Conformable format
# stopifnot(colnames(rowData(raw_merf)) == colnames(rowData(spe_visium)))

# ### Subset to Merfish/Visium overlapping hemisphere ----
# # Confirm the Pmatch column is the WM value in the turorial
# # https://jef.works/STalign/notebooks/merfish-visium-alignment-with-point-annotator.html
# hist(raw_merf$Pmatch)

# WMthresh <- 0.95
# spe_merf <- raw_merf[, raw_merf$Pmatch > WMthresh]

# Plot Visium Shape ----
spe_visium <- readRDS(
  here("processed_data", "spe_visium.rds")
)


pdf(
  here(
    "plots/case_visium_n_merf",
    "spatial_plot_raw_visium.pdf"
  ),
  height = 1 * 3, width = 1 * 3
)
spe_visium |>
  make_escheR() |>
  add_fill("sample_id", point_size = 1.5) +
  scale_fill_manual(
    values = c("visium" = "#e66101"),
    guide = "none"
  )
dev.off()
