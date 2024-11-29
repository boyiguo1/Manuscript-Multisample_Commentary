# Load data ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(escheR)
  library(here)
  library(scuttle)
  library(scran)
  library(tidyverse)
  library(sessioninfo)
})

# Load Data ----
## Visium data ----
list.files(here("raw_data/visium_mouse_brain"))

spe_visium <- read10xVisium(
  here("raw_data/visium_mouse_brain"),
  sample_id = "visium",
  type = "sparse",
  data = "filtered", # Only include in-tissue spots
  images = "lowres",
  load = FALSE
)

rowData(spe_visium)$ensmbl <- rownames(spe_visium)

### Create log counts assay ----
spe_visium <- logNormCounts(spe_visium)

## Merfish data ----
# NOTE: we used the results from Clifton et al. (2023)
merf_csv <- readr::read_csv(here("raw_data/STalign_S2R3_to_Visium.csv"))
colnames(merf_csv)[1] <- "key"
count_assay <- merf_csv[, 7:ncol(merf_csv)] |> t() # Gene-by-cell
cd <- merf_csv[, 1:6]

# Create SPE object based on data
raw_merf <- SpatialExperiment(
  assay = list("counts" = count_assay),
  sample_id = "merfish",
  colData = cd,
  spatialCoordsNames = c("STalign_x", "STalign_y")
)

rowData(raw_merf)$symbol <- rownames(raw_merf)
# Create mouse ensembl id in row data
ensmbl_vec <- rowData(spe_visium)$ensmbl[match(rowData(raw_merf)$symbol, rowData(spe_visium)$symbol)]
rowData(raw_merf)$ensmbl <- ensmbl_vec
# Non-matching ensmbl use gene symbol
rowData(raw_merf)$ensmbl[which(is.na(ensmbl_vec))] <- rowData(raw_merf)$symbol[which(is.na(ensmbl_vec))]

# Use ensmbl id as row name
rownames(raw_merf) <- rowData(raw_merf)$ensmbl

# Conformable format
stopifnot(colnames(rowData(raw_merf)) == colnames(rowData(spe_visium)))

### Subset to Merfish/Visium overlapping hemisphere ----
# Confirm the Pmatch column is the WM value in the turorial
# https://jef.works/STalign/notebooks/merfish-visium-alignment-with-point-annotator.html
hist(raw_merf$Pmatch)

WMthresh <- 0.95
spe_merf <- raw_merf[, raw_merf$Pmatch > WMthresh]

### Log counts normalizaiton ----
spe_merf <- logNormCounts(spe_merf)

spe_merf
# 40k~ cells

## Find overlapping genes ----
overlap_gene_ensmbl <- intersect(rownames(spe_visium), rownames(spe_merf))
# Error prevention
stopifnot(overlap_gene_ensmbl > 400)

## Visium with overlapping genes ----
spe_overlap_visium <- spe_visium[overlap_gene_ensmbl, ]
spe_overlap_visium <- logNormCounts(spe_overlap_visium, size.factors = NULL)

## Merfish with overlapping genes ----
spe_overlap_merf <- spe_merf[overlap_gene_ensmbl, ]
spe_overlap_merf <- logNormCounts(spe_overlap_merf, size.factors = NULL)


# Clustering and PCA (if needed) ----

## Visium (transcriptome wide) ----
set.seed(20241127)
visium_HVG_2000 <- getTopHVGs(spe_visium, n = 2000)
spe_visium <- fixedPCA(spe_visium, subset.row = visium_HVG_2000)
cluster_HVG_label <- clusterCells(spe_visium, use.dimred = "PCA")

## Visium (400 gene) ----
set.seed(20241127)
spe_overlap_visium <- fixedPCA(spe_overlap_visium, subset.row = NULL)
cluster_400_label <- clusterCells(spe_overlap_visium, use.dimred = "PCA")

## Merfish (400 gene) ----
set.seed(20241127)
spe_overlap_merf <- fixedPCA(spe_overlap_merf, subset.row = NULL)
merf_label <- clusterCells(spe_overlap_merf, use.dimred = "PCA")

## Organize cluster labels ----

spe_visium$cluster_HVG_label <- spe_overlap_visium$cluster_HVG_label <- cluster_HVG_label
spe_visium$cluster_400_label <- spe_overlap_visium$cluster_400_label <- cluster_400_label
spe_merf$merf_label <- spe_overlap_merf$merf_label <- merf_label

# plotReducedDim(spe_visium,
#   ncomponents = 4,
#   dimred = "PCA", colour_by = "cluster_HVG_label"
# )
# plotReducedDim(spe_visium,
#   ncomponents = 4,
#   dimred = "PCA", colour_by = "cluster_400_label"
# )

## Spatial Plot ----
# make_escheR(spe_visium) |>
#   add_fill(var = "cluster_HVG_label")

# make_escheR(spe_visium) |>
#   add_fill(var = "cluster_400_label")

# make_escheR(spe_overlap_merf) |>
#   add_fill(var = "merf_label", point_size = 0.5)

# Save Visium and Merfish RDS ----
saveRDS(
  spe_merf,
  here("processed_data", "spe_merf.rds")
)

saveRDS(
  spe_visium,
  here("processed_data", "spe_visium.rds")
)

saveRDS(
  spe_overlap_merf,
  here("processed_data", "spe_overlap_merf.rds")
)

saveRDS(
  spe_overlap_visium,
  here("processed_data", "spe_overlap_visium.rds")
)

# Create joint SPE with Visium and Merfish ----
colData(spe_overlap_visium) <- DataFrame(sample_id = spe_overlap_visium$sample_id)
colData(spe_overlap_merf) <- DataFrame(sample_id = spe_overlap_merf$sample_id)
spe_joint <- cbind(spe_overlap_visium, spe_overlap_merf)

## Preprocessing ----
spe_joint <- logNormCounts(spe_joint, size.factors = NULL)
set.seed(20241129)
spe_joint <- fixedPCA(spe_joint, subset.row = NULL)

## Save Joint SPE ----
saveRDS(
  spe_joint,
  here("processed_data", "spe_joint.rds")
)


# Session Info ----
sessioninfo::session_info()
