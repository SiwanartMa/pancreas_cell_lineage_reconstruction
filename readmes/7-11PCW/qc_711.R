# Perform QC 
library(tidyverse)
library(Seurat)
library(scDblFinder)
library(patchwork)
library(ggplot2)
library(SingleCellExperiment)

## Load Seurat object ----------
load("~/Desktop/researchProject/integration/originals/human7_11W.RData")
obj <- seurat

## Interacting with a Seurat object ----------

### Check dimensions
print("Dimensions of the Seurat object")
dim(obj) # 23663 41727

### Set the default assay to RNA
DefaultAssay(obj) <- "RNA"
print("Default assay:")
obj@active.assay

### Gene expression counts
print("Expression counts of the RNA assay:")
GetAssayData(obj)[1:5,1:5] # to print out expression counts

### Check the fields in the metadata
print("The columns in the metadata:")
colnames(obj[[]])

### Check colnames and rownames
print("Colnames:")
head(colnames(obj))
print("Rownames:")
head(rownames(obj))

table(obj@meta.data$orig.ident)

## Quality control of Seurat object ----------

### Doublet scoring
sce <- as.SingleCellExperiment(obj)
sce
sce <- scDblFinder(sce, samples = sce$orig.ident)
head(sce$scDblFinder.score)
table(sce$scDblFinder.class)
### Add scDblFinder outputs into Seurat object
obj <- AddMetaData(obj,
                   metadata = sce$scDblFinder.score,
                   col.name = "scDblFinder.score")
obj <- AddMetaData(obj,
                   metadata = sce$scDblFinder.class,
                   col.name = "scDblFinder.class")

### The number of singlet and doublet
table(obj[["scDblFinder.class"]])

### The number of singlet and doublet in each sample
t1 <- as.data.frame(table(obj$orig.ident, obj$scDblFinder.class))
colnames(t1) <- c("Sample", "Class", "Count")
t1

### Visualise the number of singlet and doublet in each sample
ggplot(t1, aes(x = Sample, y = Count, fill = Class)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Singlet vs Doublet Cells per Sample from W7 to W11",
       x = "Sample", y = "Total Cell Count") +
  theme_minimal()
#ggsave("~/Desktop/researchProject/integration/outputs/QC/711_scDbl.png")

# plot distribution of scDblFinder.score
VlnPlot(obj, "scDblFinder.score", group.by = "scDblFinder.class")
#ggsave("~/Desktop/researchProject/integration/outputs/QC/711scDbl_score.png")


## Quantify mitochondrial percentage ----------
# To select out mitochondrial genes
str_subset(rownames(obj), pattern = "^MT-")

# To quantify percentage of mito genes
obj[["percent_mito"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# Violin plot of percent mitochondria
VlnPlot(obj, "percent_mito", group.by = "orig.ident", pt.size = 0) +
  geom_boxplot(width=0.1, fill="white", outlier.size = 1)

table(obj$orig.ident, obj$days)

# Make violin plot
### Manually order the samples by their collected week
ordered_samples <- c("S129", 
                     "S9", "S66", "S86", "S120", 
                     "S27", "S85",
                     "S14", "S35",
                     "S13")
obj@meta.data$orig.ident <- factor(obj@meta.data$orig.ident, levels = ordered_samples)

ggplot(obj@meta.data, aes(x = orig.ident, y = percent_mito, fill = days)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.size = 1) +
  labs(title = "Percentage of mitochondrial genes", x = "Sample", y = "Mitochondrial genes%") +
  scale_fill_manual(values = c("W7" = "orange", "W8" = "#377EB8", "W9" = "#4DAF4A", "W10" = "purple", "W11" = "yellow")) +  # set custom colors
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )
#ggsave("Desktop/researchProject/integration/outputs/QC/711_mitoVln.png")


# Quantify ribosomal percentage ----------
str_subset(rownames(obj), pattern = "^RP[SL]")
obj[["percent_ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")

ggplot(obj@meta.data, aes(x = orig.ident, y = percent_ribo, fill = days)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.size = 1) +
  labs(title = "Percentage of ribosomal genes", x = "Sample", y = "Ribosomal genes%") +
  scale_fill_manual(values = c("W7" = "orange", "W8" = "#377EB8", "W9" = "#4DAF4A", "W10" = "purple", "W11" = "yellow")) +  # set custom colors
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )
#ggsave("Desktop/researchProject/integration/outputs/QC/711_riboVln.png")

### Visualise the correlation of mitochondrial genes and ribosomal genes with a scatterplot
color_list <- c("orange", #W7
                "deepskyblue4", "cornflowerblue", "deepskyblue2", #W5
                "aquamarine2", "aquamarine4","darkolivegreen4" ) #W6

FeatureScatter(obj, "percent_mito", "percent_ribo",
               group.by = "orig.ident",
               split.by = "days")

#ggsave("Desktop/researchProject/integration/outputs/QC/46_corr_ribomito.png")

## Classify cell cycle phase ----------
cc.gene <- cc.genes.updated.2019
obj <- NormalizeData(obj)  # briefly normalize data

# actual scoring of cell phases
obj <- CellCycleScoring(obj,
                        s.features = cc.genes$s.genes,
                        g2m.features = cc.genes$g2m.genes,
                        set.ident = FALSE)

# check metadata
head(obj[[]])

# tabulate the number of cell cycle phases
table(obj[["Phase"]])

## Filter out low-quality cells
day_colors <- c("W4" = "orange", "W5" = "#377EB8", "W6" = "#4DAF4A")

features <- c("percent_mito", "percent_ribo", "nCount_RNA", "nFeature_RNA")

vln_list <- VlnPlot(obj,
                    features = features,
                    group.by = "orig.ident",
                    pt.size = 0,
                    combine = FALSE)

# Add boxplot layer to each
vln_list <- lapply(vln_list, function(p) {
  p + geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
    scale_fill_manual(values = day_colors, name = "PCW") +
    aes(fill = obj@meta.data$days)  # Set fill by day for each violin
})

# Combine plots with 2 columns
wrap_plots(vln_list, ncol = 2)

#ggsave("./Desktop/researchProject/integration/outputs/QC/46_allvln.png")


## Subset the Seurat object ----------

## Subset with the criteria on paper
obj_sub<- subset(obj, percent_mito < 15 &
                   nCount_RNA > 4000 & nCount_RNA < 50000 &
                   nFeature_RNA > 500 & nFeature_RNA < 8000 &
                 scDblFinder.class == "singlet")

## Remove processed data generated by the publisher
cols_to_remove <- c("percent.mt", "old.ident", "nCount_SCT", "nFeature_SCT", "SCT_snn_res.0.5", "seurat_clusters")
for (col in cols_to_remove){
  obj_sub[[col]] <- NULL
}

## Remove the SCT assay from publisher
obj_sub[["SCT"]] <- NULL

#saveRDS(obj_sub, "~/Desktop/researchProject/integration/outputs/QC_seurat/human4_6W.rds")
