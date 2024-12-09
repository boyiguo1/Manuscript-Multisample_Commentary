# Load Library ----
library(here)
library(SpatialExperiment)
library(scater)
library(scran)
library(scuttle)
library(sessioninfo)

# Load Data ---
raw_visium <- readRDS(here("processed_data", "visium_spe.rds"))
raw_merf <- readRDS(here("processed_data", "merf_spe.rds"))

# Create a joint spe project ----
## Find overlapping genes ----
overlap_gene_symbol <- intersect(rowData(raw_visium)$symbol, rownames(raw_merf))
overlap_gene_ensmbl <- rownames(raw_visium)[match(overlap_gene_symbol, rowData(raw_visium)$symbol)]

spe_visium <- raw_visium[overlap_gene_ensmbl, ]
spe_merf <- raw_merf[overlap_gene_symbol, ]
colData(spe_visium) <- DataFrame(sample_id = spe_visium$sample_id)
colData(spe_merf) <- DataFrame(sample_id = spe_merf$sample_id)


rowData(spe_merf)$symbol <- overlap_gene_symbol
rownames(rowData(spe_merf)) <- overlap_gene_ensmbl
rownames(spe_merf) <- overlap_gene_ensmbl


spe_joint <- cbind(
  spe_visium,
  spe_merf
)

# Log normalizations ----
assay(spe_joint, "joint_logcounts") <- scuttle::normalizeCounts(spe_joint)



# Analysis ----

## Calculate PCA with individual lognorm ----
set.seed(100) 
spe_joint <- fixedPCA(spe_joint, subset.row=NULL) # Use the 450~ overlapping genes
plotReducedDim(spe_joint[, sample(ncol(spe_joint))], dimred="PCA", colour_by = "sample_id")

# TODO: Currnetly overplotting - would greatly benefitted from making each data points smaller and more transparent

plotReducedDim(spe_joint[, sample(ncol(spe_joint))], ncomponents = 4, dimred="PCA", colour_by = "sample_id")

## Calculate PCA using joint lognorm ----
set.seed(100) 
spe_joint <- fixedPCA(spe_joint, subset.row=NULL, assay.type = "joint_logcounts", name = "PCA_joint") # Use the 
plotReducedDim(spe_joint[, sample(ncol(spe_joint))], 
dimred="PCA_joint", colour_by = "sample_id")


## Clustering with joint lognorm ----
library(scran)
nn.clusters <- clusterCells(spe_joint, use.dimred="PCA_joint")
colLabels(spe_joint) <- nn.clusters

# Sesssion Info ----
sessioninfo::session_info()
