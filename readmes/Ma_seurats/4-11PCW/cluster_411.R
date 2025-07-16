# Part2: Seurat clustering

library(tidyverse)
library(Seurat)
library(clustree)

options(future.globals.maxSize= 8 * 1024^3)
mem.maxVSize(vsize = 120000)

# Load Seurat object
obj_46 <- readRDS("~/Desktop/researchProject/integration/outputs/QC_seurat/4_11PCW/human4_6W.rds")
obj_711 <- readRDS("~/Desktop/researchProject/integration/outputs/QC_seurat/4_11PCW/human7_11W.rds")

# Plot a scatter of nCount_RNA vs TUBB to see if they are positive correlate
# If yes, regress out nCount_RNA
FeatureScatter(obj_46, "nCount_RNA", "TUBB", slot = "counts")
FeatureScatter(obj_711, "nCount_RNA", "TUBB", slot = "counts")

# run SCTransform workflow
obj_46 <- SCTransform(obj_46, method = "glmGamPoi", 
                   vars.to.regress = c("nCount_RNA", "S.Score", "G2M.Score"))

obj_711 <- SCTransform(obj_711, method = "glmGamPoi", 
                      vars.to.regress = c("nCount_RNA", "S.Score", "G2M.Score"))

# Save the SCTransformed object
#saveRDS(obj_46, "~/Desktop/researchProject/integration/outputs/QC_seurat/4_11PCW/human4_6W_SCT.rds")
#saveRDS(obj_711, "~/Desktop/researchProject/integration/outputs/QC_seurat/4_11PCW/human7_11W_SCT.rds")

# Merge human7_11 with human4_6
obj <- merge(x=obj_46, y=list(obj_711))
colnames(obj[[]])
#saveRDS(obj, "~/Desktop/researchProject/integration/outputs/QC_seurat/4_11PCW/human4_11W.rds")
obj[["RNA"]]

# Find variable features
obj <- FindVariableFeatures(obj)

# PCA
obj <- RunPCA(obj)

# visualise dim reduction
DimPlot(obj)
DimPlot(obj, group.by = "orig.ident")
DimPlot(obj, group.by = "days")

# Colourcode cells by gene expression
FeaturePlot(obj, "TUBB")
FeaturePlot(obj, "percent_mito")
FeaturePlot(obj, "percent_ribo")
FeaturePlot(obj, "S.Score")
FeaturePlot(obj, "G2M.Score")

## Non-linear dimensional reduction
obj <- RunUMAP(obj, dims = 1:50)
obj <- RunTSNE(obj, dims = 1:50)

d1 <- DimPlot(obj, reduction = "umap", label = T) + NoLegend()
d2 <- DimPlot(obj, group.by = "orig.ident", reduction = "umap")
obj@meta.data$days <- factor(obj@meta.data$days, levels = c("W4", "W5", "W6", "W7", "W8", "W9", "W10", "W11"))
d3 <- DimPlot(obj, group.by = "days", reduction = "umap")
d1+d2+d3

DimPlot(obj, group.by = "orig.ident", reduction = "tsne")

### Clustering
obj <- FindNeighbors(obj, dims = 1:50)
res.range <- seq(0.5, 1.2, 0.1)
res.range

head(obj[[]])
for (res in res.range){
  obj <- FindClusters(obj, resolution = res)
}

head(obj[[]])

# Visualise clusters using clustree
clustree <- clustree(obj, prefix = "SCT_snn_res.")
clustree

# cluster using res.1
obj <- FindClusters(obj, resolution= 0.7)

DimPlot(obj, reduction = "umap", label = TRUE)

saveRDS(obj, "~/Desktop/researchProject/integration/outputs/QC_seurat/4_11PCW/human4_11W.rds")
