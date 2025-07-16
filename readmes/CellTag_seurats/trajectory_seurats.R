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
#obj <- readRDS("~/Desktop/researchProject/integration/outputs/ScaleRNA/QC_seurat/combined.rds")
obj.filt <- readRDS("~/Desktop/researchProject/integration/outputs/ScaleRNA/QC_seurat/combined_filt.rds")

# filtering out cells
setwd("~/Desktop/researchProject/integration/originals/ScaleRNA/")

# How many cell are there in each sample in the untransduced and transduced group
table(obj.filt$group, obj.filt$day)

obj.filt <- subset(obj,
                   seurat_clusters %in% c(3, 4, 7, 11))

DimPlot(obj.filt, group.by = "group", reduction = "umap", label = F)

markers <- list(
  hPSCs = c("SOX2", "NANOG"),
  Definitive_endoderm = c("FOXA2", "SOX17", "CXCR4"),
  Foregut = c("FOXA2", "SOX17", "PDX1", "HNF1B","GATA4"),
  Pancreatic_progenitor = c("PDX1", "NKX6-1", "SOX9")
)

FeaturePlot(obj.filt, unique(unlist(markers)), label = F, ncol = 5) 
p <- FeaturePlot(obj.filt,  unique(unlist(markers)), combine = FALSE)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

cowplot::plot_grid(plotlist = p, ncol = 5)

# Make Dotplot
DotPlot(obj.filt, unique(unlist(markers)), dot.scale = 5) + 
  scale_color_gradient (low = "darkslateblue", high = "brown2")+ 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 10, face = "bold"))

DotPlot(obj.filt, markers$Pancreatic_progenitor, dot.scale = 5) + 
  scale_color_gradient (low = "darkslateblue", high = "brown2")+ 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 10, face = "bold"))

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
