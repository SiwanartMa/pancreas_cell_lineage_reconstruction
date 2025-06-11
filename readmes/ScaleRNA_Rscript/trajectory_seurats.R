# Trajectory analysis

# load dependencies
library(tidyverse)
library(Seurat)
library(slingshot)
library(destiny)
library(tradeSeq)
options(future.globals.maxSize= 8 * 1024^3)
mem.maxVSize(vsize = 80000)

# Load Seurat object
obj <- readRDS("~/Desktop/researchProject/integration/outputs/ScaleRNA/QC_seurat/combined.rds")

# filtering out cells
setwd("~/Desktop/researchProject/integration/originals/ScaleRNA/")

obj.filt <- subset(obj,
                   seurat_clusters %in% c(3, 4, 7, 11))

DimPlot(obj.filt, reduction = "umap", label = T)

# Annotate clusters
new.cluster.ids <- c("Pancreatic_progenitor", "Foregut", "Definitive_endoderm", "hPSC")


names(new.cluster.ids) <- levels(obj.filt)
obj.filt <- RenameIdents(obj.filt, new.cluster.ids)
DimPlot(obj.filt, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# Pseudotime
data <- as.matrix(t(obj.filt@assays$SCT@data))
data <- data[, VariableFeatures(obj.filt)]
dm <- DiffusionMap(data, n_pcs = 3)

obj.filt[["diffmap"]] <- CreateDimReducObject(embeddings = dm@eigenvectors[,1:3], key = "dc_", assay = "SCT")
DimPlot(obj.filt, reduction = "diffmap", label = F) 

sds <- slingshot(Embeddings(obj.filt, "diffmap"), 
                 clusterLabels = obj.filt$seurat_clusters,
                 start.clus= '11')

slingLineages(sds)
slingPseudotime(sds)[,1]

obj.filt$pseudotime <- slingPseudotime(sds)[,1]

FeaturePlot(obj.filt, "pseudotime", reduction = "umap", label = T)

saveRDS(obj.filt, "~/Desktop/researchProject/integration/outputs/ScaleRNA/QC_seurat/combined_filt.rds")
