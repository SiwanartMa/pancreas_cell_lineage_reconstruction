library(Seurat)
library(ggplot2)
library(patchwork)
library(SeuratWrappers)
mem.maxVSize(vsize = 80000)

obj <- readRDS("./Desktop/researchProject/integration/outputs/integrated_interactive.rds")
obj

obj <- FindNeighbors(obj, reduction = "integrated.pca", dims = 1:30)

obj <- FindClusters(obj, resolution = 2, cluster.name = "cca_clusters")

obj <- RunUMAP(obj, reduction = "integrated.pca", dims = 1:30, reduction.name = "umap.cca")

p1_days <- DimPlot(
  obj,
  reduction = "umap.cca",
  group.by = c("days"),
  combine = FALSE, label.size = 2
)
p1_days
ggsave(file="./Desktop/researchProject/integration/outputs/integrated_figure/umap_days.png",
       width=8,
       height=6,
       dpi=300)

p1_time <- DimPlot(
  obj,
  reduction = "umap.cca",
  group.by = c("time"),
  combine = FALSE, label.size = 2
)
p1_time
ggsave(file="./Desktop/researchProject/integration/outputs/integrated_figure/umap_time.png",
       width=8,
       height=6,
       dpi=300)


p1 <- DimPlot(
  obj,
  reduction = "umap.cca",
  group.by = "cca_clusters",
  combine = FALSE, label.size = 2
)
p1
ggsave(file="./Desktop/researchProject/integration/outputs/integrated_figure/umap_ccaclusters.png",
       width=8,
       height=6,
       dpi=300)
