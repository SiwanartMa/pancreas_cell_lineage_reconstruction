# Trajectory analysis

# load dependencies
library(tidyverse)
library(Seurat)
library(slingshot)
library(destiny)
library(tradeSeq)
library(harmony)
library(clustree)
library(EnhancedVolcano)
library(patchwork)
options(future.globals.maxSize= 8 * 1024^3)
mem.maxVSize(vsize = 120000)

# Import the seurat object
obj.filt <- readRDS("~/Desktop/researchProject/integration/outputs/QC_seurat/4_11PCW/epith_411.rds")

## create a diffusion map
data <- GetAssayData(obj.filt, assay = "SCT", layer = "data")

data <- as.matrix(t(data))

data <- data[, VariableFeatures(obj.filt)]
dm <- DiffusionMap(data, n_pcs = 3)
obj.filt[["diffmap"]] <- CreateDimReducObject(embeddings = dm@eigenvectors,
                                              key = "dc_",
                                              assay = "SCT")

DimPlot(obj.filt, reduction = "diffmap", label = T)

# run slingshot
# start from dorsal MP
sds <- slingshot(Embeddings(obj.filt, "diffmap"),
                 clusterLabels = obj.filt$seurat_clusters, 
                 start.clus = "12")

slingLineages(sds)
slingPseudotime(sds)

# start from PB
sds <- slingshot(Embeddings(obj.filt, "diffmap"),
                 clusterLabels = obj.filt$seurat_clusters, 
                 start.clus = "16")

slingLineages(sds)
slingPseudotime(sds)[,1]

obj.filt$pseudotime <- slingPseudotime(sds)[,1]

FeaturePlot(obj.filt, "pseudotime", reduction = "umap", label = T,) +
  scale_color_gradient(low = "lightgreen", high = "orange")

# Not save yet
#vsaveRDS(obj.filt,"~/Desktop/researchProject/integration/outputs/QC_seurat/4_11PCW/epith_411.rds" )

