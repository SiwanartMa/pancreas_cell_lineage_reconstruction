library(Seurat)
library(glue)
library(patchwork)
library(BPCells)
library(scDblFinder)
library(tidyverse)
library(dplyr)

# import the celltagged Seurat objects have been processed with scDblFinder
for (i in c(2, 6, 8, 10)) {
  path <- glue("/Users/mayongzhi/Desktop/researchProject/integration/outputs/ScaleRNA/QC_seurat/d{i}_tag.rds")
  obj <- readRDS(path)
  assign(glue("d{i}_tag"), obj)
}

for (i in c(2, 6, 8, 10)) {
  path <- glue("/Users/mayongzhi/Desktop/researchProject/integration/outputs/ScaleRNA/QC_seurat/d{i}.rds")
  obj <- readRDS(path)
  assign(glue("d{i}"), obj)
}

ctrl_tag <- readRDS("/Users/mayongzhi/Desktop/researchProject/integration/outputs/ScaleRNA/QC_seurat/ctrl_tag.rds")
ctrl <- readRDS("/Users/mayongzhi/Desktop/researchProject/integration/outputs/ScaleRNA/QC_seurat/ctrl.rds")

# Merge all objects
obj <- merge(x = ctrl_tag, y = list(d2_tag, d6_tag, d8_tag, d10_tag,
                                    ctrl, d2, d6, d8, d10))
colnames(obj[[]])

# Rename samples by removing ".ScaleRNA" from the sample column
obj@meta.data$sample <- gsub("\\.ScaleRNA", "", obj@meta.data$sample)


# To quantify percentage of mito genes
setwd("~/Desktop/researchProject/integration/originals/ScaleRNA/")
obj[["percent_mito"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# Violin plot of percent mitochondria
ordered_samples <- c("iPSCs-Frozen", "CTRL-D2", "CTRL-D6", "CTRL-D8", "CTRL-D10",
                     "CellTag-iPSCs-Frozen", "CellTag-D2", "CellTag-D6", "CellTag-D8", "CellTag-D10")
                     
obj@meta.data$sample <- factor(obj@meta.data$sample, levels = ordered_samples)

VlnPlot(obj, "percent_mito", group.by = "sample", pt.size = 0) +
  geom_boxplot(width=0.1, fill="white", outlier.size = 1)

# To quantify percentage of ribo genes
str_subset(rownames(obj), pattern = "^RP[SL]")
obj[["percent_ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")

VlnPlot(obj, "percent_ribo", group.by = "sample", pt.size = 0) +
  geom_boxplot(width=0.1, fill="white", outlier.size = 1)

# Plot the correlation between percent_mito and percent_ribo
FeatureScatter(obj, "percent_mito", "percent_ribo",
               group.by = "sample")

# Identify cell phase
cc.gene <- cc.genes.updated.2019
obj <- NormalizeData(obj) 

obj <- CellCycleScoring(obj,
                        s.features = cc.genes$s.genes,
                        g2m.features = cc.genes$g2m.genes,
                        set.ident = FALSE)
table(obj[["Phase"]])

# Plot all QC features
features <- c("nCount_RNA", "nFeature_RNA")

vln_list <- VlnPlot(obj,
                    features = features,
                    group.by = "sample",
                    pt.size = 0, ) + NoLegend()

# Add boxplot layer to each
vln_list <- lapply(vln_list, function(p) {
  p + geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5)
})

# Combine plots with 2 columns
wrap_plots(vln_list, ncol = 2)

# Subset cells
obj <- subset(obj, 
              percent_mito < 5 &
              nFeature_RNA > 500 &
              scDblFinder.class == "singlet")

# Save the seurat object
#saveRDS(obj, "~/Desktop/researchProject/integration/outputs/ScaleRNA/QC_seurat/combined.rds")
