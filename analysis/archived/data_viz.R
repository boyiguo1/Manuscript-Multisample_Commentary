plotReducedDim(spe_visium,
  ncomponents = 4,
  dimred = "PCA", colour_by = "cluster_HVG_label"
)
plotReducedDim(spe_visium,
  ncomponents = 4,
  dimred = "PCA", colour_by = "cluster_400_label"
)


## Spatial Plot ----
make_escheR(spe_visium) |>
  add_fill(var = "cluster_HVG_label")

make_escheR(spe_visium) |>
  add_fill(var = "cluster_400_label")

make_escheR(spe_overlap_merf) |>
  add_fill(var = "merf_label", point_size = 0.5)


# Descriptive Analysis ----
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

## Lib size of overlapping genes ----
overlap_gene_symbol <- intersect(
  rowData(spe_visium)$symbol, rownames(sub_merf)
)
overlap_gene_ensmbl <- rownames(spe_visium)[
  match(overlap_gene_symbol, rowData(spe_visium)$symbol)
]
overlap_merf_lib_size <- counts(sub_merf)[overlap_gene_symbol, ] |>
  colSums()
overlap_visium_lib_size <- counts(spe_visium)[overlap_gene_ensmbl, ] |>
  colSums()
rbind(
  data.frame(
    lib_size = merf_lib_size
  ) |> mutate(tech = "merfish"),
  data.frame(
    lib_size = visium_lib_size
  ) |> mutate(tech = "visium"),
  data.frame(
    lib_size = overlap_merf_lib_size
  ) |> mutate(tech = "merfish_overlap"),
  data.frame(
    lib_size = overlap_visium_lib_size
  ) |> mutate(tech = "visium_overlap")
) |>
  ggplot() +
  stat_ecdf(
    aes(x = lib_size, color = tech),
    geom = "step"
  ) +
  scale_x_log10() +
  theme_minimal()

# TODO:
# 1) make the x_axis more readable.

## gene ----
gene_symbol <- "Baiap2"
em_id <- rownames(rowData(spe_visium))[which(rowData(spe_visium)$symbol == gene_symbol)]
lc_merf <- logcounts(sub_merf)[gene_symbol, ]
lc_visium <- logcounts(spe_visium)[em_id, ]

lc_visium_overlap <- normalizeCounts(spe_visium,
  size.factors = NULL,
  subset.row = overlap_gene_ensmbl
)[em_id, ]
lc_merf_overlap <- normalizeCounts(sub_merf,
  size.factors = NULL, subset.row = overlap_gene_symbol
)[gene_symbol, ]

# calculateSumFactors()

rbind(
  data.frame(
    lc = lc_merf
  ) |> mutate(tech = "merfish"),
  data.frame(
    lc = lc_visium
  ) |> mutate(tech = "visium"),
  data.frame(
    lc = lc_merf_overlap
  ) |> mutate(tech = "merfish_overlap"),
  data.frame(
    lc = lc_visium_overlap
  ) |> mutate(tech = "visium_overlap")
) |>
  ggplot() +
  stat_ecdf(
    aes(x = lc, color = tech),
    geom = "step"
  ) +
  # scale_x_log10() +
  theme_minimal() +
  labs(title = gene_symbol)

# Create joint data
all_spe <- c(visium_spe, sub_merf)

# Plots ----
# make_escheR(spe_merf) |>
#   add_fill("sample_id", point_size = 0.2)

# make_escheR(
#   data.frame(sample_id = "merfish"),
#   .x = merf_csv$center_x, .y = merf_csv$center_y
# ) |>
#   add_fill("sample_id", point_size = 0.2)
# make_escheR(
#   data.frame(sample_id = "merfish"),
#   .x = sub_merf$center_x, .y = sub_merf$center_y
# ) |>
#   add_fill("sample_id", point_size = 0.2)
