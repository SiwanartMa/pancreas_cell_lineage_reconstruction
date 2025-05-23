library(Seurat)
library(glue)
library(patchwork)
library(BPCells)
library(scDblFinder)
library(tidyverse)
#install.packages('BPCells', repos = c('https://bnprks.r-universe.dev', 'https://cloud.r-project.org'))


# import the celltagged Seurat objects
for (i in c(2, 6, 8, 10)) {
  path <- glue("/Users/mayongzhi/Desktop/researchProject/integration/originals/ScaleRNA/CellTag-D{i}.ScaleRNA_SeuratObject.rds")
  obj <- readRDS(path)
  assign(glue("d{i}_tag"), obj)
}

ctrl_tag <- readRDS("/Users/mayongzhi/Desktop/researchProject/integration/originals/ScaleRNA/CellTag-iPSCs-Frozen.ScaleRNA_SeuratObject.rds")

obj <- merge(x = ctrl_tag, y = list(d2_tag, d6_tag, d8_tag, d10_tag))

# To quantify percentage of mito genes
obj[["percent_mito"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# Violin plot of percent mitochondria
ordered_samples <- c("CellTag-iPSCs-Frozen.ScaleRNA", "CellTag-D2.ScaleRNA", "CellTag-D6.ScaleRNA", "CellTag-D8.ScaleRNA", "CellTag-D10.ScaleRNA")
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
features <- c("percent_mito", "percent_ribo", "nCount_RNA", "nFeature_RNA")

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

# Get count matrix
raw_counts <- as.matrix(GetAssayData(obj, layer = "counts"))
sce <- scDblFinder(, samples = obj$sample)

obj <- AddMetaData(obj,
                   metadata = sce$scDblFinder.score,
                   col.name = "scDblFinder.score")
obj <- AddMetaData(obj,
                   metadata = sce$scDblFinder.class,
                   col.name = "scDblFinder.class")

table(obj[["scDblFinder.class"]])

VlnPlot(d2_tag, "scDblFinder.score", group.by = "scDblFinder.class")

table(obj$RT)
table(obj$i5)
table(obj$Ligation)
