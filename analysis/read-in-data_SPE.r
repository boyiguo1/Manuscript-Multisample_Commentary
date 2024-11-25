# Load data ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(escheR)
  # library(SpatialFeatureExperiment)
  library(here)
  library(sessioninfo)
  library(scuttle)
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

spe_visium
# 2264 spots and 19465 genes

### Create log counts assay ----

spe_visium <- logNormCounts(spe_visium)

## Merfish data ----
# sfe_mer <- readVizgen(here("raw_data/merfish_slice2_rep3"),
# min_area = 0,  z = "all")


### Subset to overlapping regions only ----
# NOTE: we used the results from Clifton et al. (2023)
merf_csv <- readr::read_csv(here("raw_data/STalign_S2R3_to_Visium.csv"))

colnames(merf_csv)[1] <- "key"


count_assay <- merf_csv[, 7:ncol(merf_csv)] |> t() # Gene-by-cell

cd <- merf_csv[, 1:6]

spe_merf <- SpatialExperiment(
  assay = list("counts" = count_assay),
  colData = cd,
  spatialCoordsNames = c("STalign_x", "STalign_y")
)
spe_merf$sample_id <- "merfish"

# Confirm the Pmatch column is the WM value in the turorial
# https://jef.works/STalign/notebooks/merfish-visium-alignment-with-point-annotator.html
hist(spe_merf$Pmatch)

make_escheR(spe_merf) |>
  add_fill("sample_id", point_size = 0.2)

make_escheR(
  data.frame(sample_id = "merfish"),
  .x = merf_csv$center_x, .y = merf_csv$center_y
) |>
  add_fill("sample_id", point_size = 0.2)


WMthresh <- 0.95
sub_merf <- spe_merf[, spe_merf$Pmatch > WMthresh]

make_escheR(
  data.frame(sample_id = "merfish"),
  .x = sub_merf$center_x, .y = sub_merf$center_y
) |>
  add_fill("sample_id", point_size = 0.2)

sub_merf
# 40k~ cells


### Log counts normalizaiton ----
sub_merf <- logNormCounts(sub_merf)


# Save RDS ----




# Analysis ----
## Overlapping genes ----
intersect(rowData(spe_visium)$symbol, rownames(sub_merf))

# TODO
# 1) do we need visualization for this?
# 2) if yes, what forms to highlight what knowledg
#    * Ven Diaggram to highlight different sets of genes
#    * Scatter plot for two different Resolution trade-off
#    * Schematic

### Venn Diagram ----
library(VennDiagram)

dev.off()
venn.diagram(
  x = list(
    rowData(spe_visium)$symbol,
    rownames(sub_merf)
  ),
  filename = NULL,
  category.names = c("Visium", "Merfish"),
  output = FALSE
) |> grid.draw()

# TODO: aethestic adjustment
# https://r-graph-gallery.com/14-venn-diagramm


### Scatter plot ----
p_df <- data.frame(
  tech = c("Visium", "Merfish"),
  gene = c(nrow(spe_visium), nrow(sub_merf)),
  sample_size = c(ncol(spe_visium), ncol(sub_merf))
)

plot(p_df$gene, p_df$sample_size)
# TODO:
# 1) log scaling for x and y
# 2) add label to the data points

# length(rownames(sub_merf))

## Library size ----
merf_lib_size <- counts(sub_merf) |> colSums()
visium_lib_size <- counts(spe_visium) |> colSums()

library(tidyverse)

rbind(
  data.frame(
    lib_size = merf_lib_size
  ) |> mutate(tech = "merfish"),
  data.frame(
    lib_size = visium_lib_size
  ) |> mutate(tech = "visium")
) |>
ggplot() +
  stat_ecdf(aes(x = lib_size, color = tech), geom = "step") +
  scale_x_log10() +
  theme_minimal()

# TODO: 
# 1) make the x_axis more readable.

## gene ----
gene_symbol <- "Baiap2"
em_id <- which(rowData(spe_visium)$symbol == gene_symbol)
lc_merf <- logcounts(sub_merf)[gene_symbol, ]
lc_visium <- logcounts(spe_visium)[em_id, ]

rbind(
  data.frame(
    lc = lc_merf
  ) |> mutate(tech = "merfish"),
  data.frame(
    lc = lc_visium
  ) |> mutate(tech = "visium")
) |>
ggplot() +
  stat_ecdf(aes(x = lc, color = tech), geom = "step") +
  # scale_x_log10() +
  theme_minimal() +
  labs(title = gene_symbol)


## Dimension reduction


# Session Info ----
sessioninfo::session_info()
