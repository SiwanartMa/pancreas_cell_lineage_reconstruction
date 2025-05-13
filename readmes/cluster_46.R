# Part2: Seurat clustering

library(tidyverse)
library(Seurat)
library(clustree)

options(future.globals.maxSize= 8 * 1024^3)
mem.maxVSize(vsize = 80000)
# Load Seurat object
obj <- readRDS("~/Desktop/researchProject/integration/outputs/QC_seurat/human4_6W.rds")

# plot a scatter of nCount_RNA vs TUBB
FeatureScatter(obj, "nCount_RNA", "TUBB", slot = "counts")
# It is positive correlate with the sequencing depth, which is bad

# run SCTransform workflow
obj <- SCTransform(obj, method = "glmGamPoi")

# check again scatterplot
FeatureScatter(obj, "nCount_RNA", "TUBB", slot = "counts")
# Not positive correlate anymore, which is good!

### Linear Dimension Reduction

# PCA
obj <- RunPCA(obj)
# extract out the best 50 dimensions and save it in a new slot

# visualise dim reduction
DimPlot(obj)
DimPlot(obj, group.by = "orig.ident")

# Colourcode cells by gene expression
FeaturePlot(obj, "TUBB")
FeaturePlot(obj, "percent_mito")
FeaturePlot(obj, "percent_ribo")
FeaturePlot(obj, "S.Score")
FeaturePlot(obj, "G2M.Score")

## Non-linear dimensional reduction
obj <- RunUMAP(obj, dims = 1:50)
obj <- RunTSNE(obj, dims = 1:50)


DimPlot(obj, group.by = "orig.ident", reduction = "umap")
DimPlot(obj, group.by = "orig.ident", reduction = "tsne")

### Clustering
obj <- FindNeighbors(obj, dims = 1:50)

res.range <- seq(0.1, 1,0.1)
res.range

head(obj[[]])
for (res in res.range){
  obj <- FindClusters(obj, resolution = res)
}

head(obj[[]]) # The column seurat_clusters will be the last res (2) due to the for loop

# Visualise clusters using clustree
clustree <- clustree(obj, prefix = "SCT_snn_res.") # choose
ggsave("~/Desktop/researchProject/integration/outputs/QC/46_clustree.png")

# cluster using res.1.6
obj <- FindClusters(obj, resolution=1)

DimPlot(obj, reduction = "umap", label = TRUE)

saveRDS(obj, "~/Desktop/researchProject/integration/outputs/QC_seurat/human4_6W.rds")
