# NOTE: for some reason, SFE can not load the MEFISH data correctly, as it only load in 2k cells.
# Deprecated code for now. 
# Load data ----
suppressPackageStartupMessages({
  library(SpatialFeatureExperiment)
  library(here)
  library(sessioninfo)
})

# Load Data ----
## Visium data ----
list.files(here("raw_data/visium_mouse_brain"))

sfe_visium <- read10xVisiumSFE(
  here("raw_data/visium_mouse_brain"),
  sample_id = "visium",
  type = "sparse",
  data = "filtered",
  images = "lowres",
  zero.policy = TRUE,
)

sfe_visium
#2264 spots and 19465 genes

## Merfish data ----
sfe_mer <- readVizgen(here("raw_data/merfish_slice2_rep3"),
min_area = 0,  z = "all")


### Subset to overlapping regions only ----
# NOTE: we used the results from Clifton et al. (2023)
subset_id <- readr::read_csv(here("raw_data/STalign_S2R3_to_Visium.csv"), col_select = c(1),
col_types = "c")

read_csv("{file_path}", col_types = cols(X1 = col_character()))

# Save RDS ----



# Session Info ----
sessioninfo::session_info()
