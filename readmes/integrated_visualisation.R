library(Seurat)
library(ggplot2)
library(patchwork)
#library(SeuratWrappers)
mem.maxVSize(vsize = 80000)

# Load the integrated Seurat object
obj <- readRDS("/Users/mayongzhi/Desktop/researchProject/integration/outputs/integrated_interactive.rds")

# The dimension and resolution is based on Ma et al. 2023
obj <- FindNeighbors(obj, reduction = "integrated.pca", dims = 1:50)
obj <- FindClusters(obj, resolution = 1, cluster.name = "cca_clusters")
obj <- RunUMAP(obj, reduction = "integrated.pca", dims = 1:50, reduction.name = "umap.cca")

# Plot
p1_days <- DimPlot(
  obj,
  reduction = "umap.cca",
  group.by = c("days"),
  combine = FALSE, label.size = 2
)
p1_days

#ggsave(file="./Desktop/researchProject/integration/outputs/integrated_figure/umap_days.png",
#       width=8,
#       height=6,
#       dpi=300)

p1_time <- DimPlot(
  obj,
  reduction = "umap.cca",
  group.by = c("time"),
  combine = FALSE, label.size = 2
)
p1_time


p1 <- DimPlot(
  obj,
  reduction = "umap.cca",
  group.by = "cca_clusters",
  combine = FALSE, label.size = 2
)
p1

# Do dimensional reduction to the unintegrated object to see the influence of integration
merged_obj = readRDS("/Users/mayongzhi/Desktop/researchProject/integration/outputs/merged_object_pca.rds")
merged_obj <- FindNeighbors(merged_obj,reduction = "pca", dims = 1:50)
merged_obj <- FindClusters(merged_obj, resolution = 1, cluster.name = "pca_clusters")
merged_obj <- RunUMAP(merged_obj, reduction = "pca", dims = 1:50, reduction.name = "umap.pca")

# Plot
p1_merged_time <- DimPlot(
  merged_obj,
  reduction = "umap.pca",
  group.by = c("time"),
  combine = FALSE, label.size = 2
)
p1_time
