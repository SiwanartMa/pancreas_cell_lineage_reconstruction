# 4- Differential gene expression

# load packages
library(tidyverse)
library(Seurat)
library(patchwork)
options(future.globals.maxSize= 8 * 1024^3)
mem.maxVSize(vsize = 80000)

# Load Seurat object
obj <- readRDS("~/Desktop/researchProject/integration/outputs/QC_seurat/human4_6W.rds")

# FindAllMarkers() was performed on HPC
#diff.markers <- FindAllMarkers(obj, test.use = "wilcox", only.pos = TRUE)
diff.markers <- read.csv("~/Desktop/researchProject/integration/metadata/4_6_markers.csv")


# Generate the marker dict from paper
markers<- list(
  Epithelial = c("EPCAM", "PDX1"),
  Mesenchymal = c("COL3A1"),
  Endothelial = c("PECAM1"),
  Neural = c("ASCL1"),
  Immune = c("PTPRC"),
  Erythroid = c("HBA1"), 
  Liver = c("ALB", "AFP")
)

# Generate FeaturePlots with both cell type and gene name in the title
plot_list <- lapply(names(markers), function(celltype) {
  gene <- markers[[celltype]][1]
  FeaturePlot(obj, features = gene) +
    ggtitle(paste(celltype, "|", gene))
})

# Combine them using patchwork's `wrap_plots`
wrap_plots(plotlist = plot_list, ncol = 3)
#ggsave("~/Desktop/researchProject/integration/outputs/QC/46_celltypes.png")

DimPlot(obj, reduction="umap", group.by = "seurat_clusters", label = T) + NoLegend()
#ggsave("~/Desktop/researchProject/integration/outputs/QC_plot/4_6PCW_reQC/46_seurat_clusters.png")

FeaturePlot(obj, features = unlist(markers), label = T, ncol = 3) 
DotPlot(obj, features = markers, assay = "RNA")
DoHeatmap(obj, features = top$gene, size = 3)

# Feature plot for epithelial markers EPCAM and PDX1 
FeaturePlot(obj, c("EPCAM", "PDX1"), label = TRUE)
FeaturePlot(obj, features = c("EPCAM", "PDX1"), blend = TRUE)
#ggsave("~/Desktop/researchProject/integration/outputs/QC/46_epithelial.png")

top <- diff.markers %>% 
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = )

top$gene
FeaturePlot(obj, top$gene, label = T, ncol = 3)
#ggsave("~/Desktop/researchProject/integration/outputs/QC/46_FeaSca_topgenes.png")


DoHeatmap(obj, features = top$gene)

VlnPlot(obj, features = top$gene[1:4], ncol = 1)

# attach the diff.markers df to the Seurat object
Misc(obj, slot = "diff.markers") <- diff.markers
saveRDS(obj, "~/Desktop/researchProject/integration/outputs/QC_seurat/human4_6W.rds")
