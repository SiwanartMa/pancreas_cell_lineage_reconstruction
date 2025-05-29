#  Differential gene expression

# load packages
library(tidyverse)
library(Seurat)
library(patchwork)
library(dplyr)
options(future.globals.maxSize= 8 * 1024^3)
mem.maxVSize(vsize = 80000)

# Load Seurat object
obj <- readRDS("~/Desktop/researchProject/integration/outputs/ScaleRNA/QC_seurat/combined.rds")

# Find markers
setwd("~/Desktop/researchProject/integration/originals/ScaleRNA/")
obj <- PrepSCTFindMarkers(obj, assay = "SCT", verbose = TRUE)
#diff.markers <- FindAllMarkers(obj, test.use = "wilcox", only.pos = TRUE)
diff.markers <- read.csv("~/Desktop/researchProject/integration/metadata/scale_combined_diff_markers.csv")
#write.csv(diff.markers, "~/Desktop/researchProject/integration/metadata/scale_combined_diff_markers.csv")

markers <- list(
  hPSCs = c("SOX2", "OCT4", "NANOG"),
  Definitive_endoderm = c("FOXA2", "SOX17", "CXCR4"),
  Foregut = c("FOXA2", "SOX17", "PDX1", "HNF1B", "HNF6", "GATA4"),
  Pancreatic_progenitor = c("PDX1", "NKX6-1", "SOX9", "PTF1A"),
  Endocrine_progenitor = c("PDX1", "NKX6-1", "NGN3", "NEUROD1")
) # OCT4, HNF6, PTF1A, NGN3, NEUROD1 are not found

obj <- AddModuleScore(obj,
                      features = markers,
                      name = "Score")

obj@meta.data <- obj@meta.data %>%  
  dplyr::rename(hPSCs_Score = Score1, # hPSCs
                DE_Score = Score2, # Definitive endoderm
                F_Score = Score3, # Foregut
                PP_Score = Score4,# Pancreatic progenitor
                EP_Score = Score5) # Endocrine progenitor

# Plot feature plot split by group
features <- c("hPSCs_Score", "DE_Score", "F_Score", "PP_Score", "EP_Score")
obj$group <- factor(obj$group, levels = c("CTRL", "CellTag"))
FeaturePlot(obj, features, label = TRUE, 
            min.cutoff = 0, reduction = "umap", split.by = c("group"))

# Plot feature plot split by day
obj$day <- factor(obj$day, levels = c("iPSC", "D2", "D6", "D8", "D10"))
FeaturePlot(obj, features, label = T, 
            min.cutoff = 0, reduction = "umap", split.by = c("day"))

# Plot dot plot
DotPlot(obj, features, split.by = "group", group.by = "day", cols = c("lightgreen", "orange"))
DotPlot(obj, features, split.by = "group", cols = c("lightgreen", "orange"))

# Subset and generate the two DotPlots
d1 <- DotPlot(subset(obj, group == "CellTag"), group.by = "day", features = features) + 
  ggtitle("CellTag") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

d2 <- DotPlot(subset(obj, group == "CTRL"),  group.by = "day", features = features) + 
  ggtitle("CTRL") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine the two plots with a shared legend and overall title
combined <- d1 + d2 + 
  plot_layout(guides = "collect") & 
  plot_annotation(title = "DotPlot Comparison: CellTag vs CTRL in each stage") &
  theme(legend.position = "right")

print(combined)

# Generate individual violin plots with median dot and reference line
violin_plots <- lapply(features, function(score) {
  VlnPlot(obj, features = score, pt.size = 0) +
    stat_summary(fun = median, geom = "point", size = 2, colour = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    NoLegend() +
    ggtitle(score)
})

# Combine plots in 2 rows
combined_plot <- wrap_plots(violin_plots, ncol = 2)
combined_plot

Misc(obj, slot = "diff.markers") <- diff.markers
saveRDS(obj, "~/Desktop/researchProject/integration/outputs/ScaleRNA/QC_seurat/combined.rds")

