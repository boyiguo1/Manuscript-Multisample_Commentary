# Load library ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(here)
  library(ggplot2)
  library(ggridges)
  library(sessioninfo)
  library(escheR)
})

# Load data ----
# TODO: add this part
# Roughly needing the four spe datasets


# Library Size Plots ----
## create lib size df ----
# Calcualte lib size for each  tech
merf_lib_size <- counts(spe_merf) |>
  colSums()
visium_lib_size <- counts(spe_visium) |>
  colSums()
overlap_merf_lib_size <- counts(spe_overlap_merf) |>
  colSums()
overlap_visium_lib_size <- counts(spe_overlap_visium) |>
  colSums()

## create dataframe ----
lib_size_df <- rbind(
  data.frame(
    lib_size = merf_lib_size,
    tech = "Merfish",
    overlap = FALSE
  ),
  data.frame(
    lib_size = visium_lib_size,
    tech = "Visium",
    overlap = FALSE
  ),
  # NOTE: Decide to remove this because redundent information
  # data.frame(
  #   lib_size = overlap_merf_lib_size,
  #   tech = "merfish",
  #   overlap = TRUE
  # ),
  data.frame(
    lib_size = overlap_visium_lib_size,
    tech = "Visium\n(overlap)",
    overlap = TRUE
  )
)


## Ridge Plots ----

ridge_p <- lib_size_df |>
  ggplot(aes(x = lib_size, y = tech, fill = tech)) +
  geom_density_ridges(
    alpha = 0.5, show.legend = FALSE,
    # scale = 4, # more over lapping between ridges
    # rel_min_height = 0.001, # remove trailing tails
    # Add rug
    jittered_points = TRUE, point_shape = "|",
    position = position_points_jitter(width = 0.05, height = 0),
    point_size = 1.5
  ) +
  scale_fill_manual(
    values = c("Merfish" = "#5e3c99", "Visium\n(overlap)" = "#fdb863", "Visium" = "#e66101")
  ) +
  scale_y_discrete(
    limits = c("Visium", "Visium\n(overlap)", "Merfish"),
    expand = expand_scale(mult = c(0.01, .3))
  ) +
  scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  theme_minimal() +
  labs(
    x = "Number of Transcripts" # , y = "Technology & Overlap Status"
  ) +
  # coord_cartesian(clip = "off") +
  theme(
    axis.title.y = element_blank() # ,
    # text = element_text(size = 12)
  )



# TODO: find color for these four levels:alpha
# Arragne the order of the tech in the plot for the plot
# Doesn't have to imclude the merf_overlap
#



## Cumulative Distribution Functions ----
ecd_p <- lib_size_df |>
  ggplot() +
  stat_ecdf(
    aes(x = lib_size, color = tech),
    geom = "step", size = 2
  ) +
  scale_color_manual(
    values = c("Merfish" = "#5e3c99", "Visium\n(overlap)" = "#fdb863", "Visium" = "#e66101"),
    limits = c("Visium", "Visium\n(overlap)", "Merfish"),
  ) +
  scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  theme_minimal() +
  labs(
    x = "Number of Transcripts", y = "Empirical Cumulative Distribution"
  ) +
  theme(text = element_text(size = 12))

## Visium Spatial Plot ----
spe_visium$lib_size <- visium_lib_size
make_escheR(spe_visium) |>
  add_fill("lib_size")


## Merfish Spatial Plot ----
spe_merf$lib_size <- merf_lib_size
make_escheR(spe_merf) |>
  add_fill("lib_size", point_size = 0.5)


# Save plot ----
pdf(
  here(
    "plots/case_visium_n_merf",
    "ridge_lib_size.pdf"
  ) # ,
  # width = 1.7,
  # height = 1.7
)
# png(
#   here(
#     "plots/case_visium_n_merf",
#     "ridge_lib_size.png"
#   )#,
#   # height = 150,
#   # width = 150
#   # # units = "in"
# )
print(p)
dev.off()


# Session Info ----
session_info()
