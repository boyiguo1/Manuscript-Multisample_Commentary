# Load library ----
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(here)
  library(ggplot2)
  library(ggridges)
  library(sessioninfo)
  library(escheR)
  library(scale)
  library(viridis)
})

# Load data ----
spe_merf <- readRDS(
  here("processed_data", "spe_merf.rds")
)

spe_visium <- readRDS(
  here("processed_data", "spe_visium.rds")
)

spe_overlap_merf <- readRDS(
  here("processed_data", "spe_overlap_merf.rds")
)

spe_overlap_visium <- readRDS(
  here("processed_data", "spe_overlap_visium.rds")
)


# Library Size Plots ----
gene_symbol <- "Baiap2"
em_id <- rownames(rowData(spe_visium))[
  which(rowData(spe_visium)$symbol == gene_symbol)
]

## create dataframe ----
lib_size_df <- rbind(
  data.frame(
    lib_size = logcounts(spe_merf)[em_id, ],
    tech = "Merfish",
    overlap = FALSE
  ),
  data.frame(
    lib_size = logcounts(spe_visium)[em_id, ],
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
    lib_size = logcounts(spe_overlap_visium)[em_id, ],
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
  # scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Log2 Norm-Counts" # , y = "Technology & Overlap Status"
  ) +
  # coord_cartesian(clip = "off") +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
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
  # scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  theme_minimal() +
  labs(
    x = "Log2 Norm-Counts", y = "Empirical Cumulative Distribution"
  ) +
  theme(text = element_text(size = 12))

## Visium Spatial Plot ----
pdf(
  here(
    "plots/case_visium_n_merf",
    "spatial_plot_Baiap2_visium.pdf"
  ),
  width = 4.8,
  height = 4.2
)
spe_visium$lib_size <- logcounts(spe_visium)[em_id, ]
make_escheR(spe_visium) |>
  add_fill("lib_size", point_size = 2.2) +
  # theme_blank()+
  theme_void(base_size = 12) +
  scale_fill_viridis_c(
    # trans = "log",
    option = "inferno",
    name = "" # ,
    # breaks = trans_breaks("log10", function(x) 10^x),
    # labels = trans_format("log10", math_format(10^.x))
  )
dev.off()


## Merfish Spatial Plot ----
pdf(
  here(
    "plots/case_visium_n_merf",
    "spatial_plot_Baiap2_merfish.pdf"
  ),
  width = 4.8,
  height = 4.2
)
spe_merf$lib_size <- logcounts(spe_merf)[em_id, ]
make_escheR(spe_merf) |>
  add_fill("lib_size", point_size = 0.5) +
  # theme_blank()+
  theme_void(base_size = 12) +
  scale_fill_viridis_c(
    # trans = "log",
    option = "inferno",
    name = "" # ,
    # breaks = trans_breaks("log10", function(x) 10^x),
    # labels = trans_format("log10", math_format(10^.x))
  )
dev.off()

# Save plot ----
pdf(
  here(
    "plots/case_visium_n_merf",
    "ridge_plot_Baiap2.pdf"
  ),
  width = 3.3,
  height = 2
)
print(ridge_p)
dev.off()


# Session Info ----
session_info()
