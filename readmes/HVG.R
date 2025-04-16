# Install packages
library(Seurat)
library(tidyverse)
library(patchwork)

# Load the 3 Seurat objects
# Increase the memory
mem.maxVSize(vsize = 60000)

# Load 4-6 weeks Seurat object
load("/Users/mayongzhi/Desktop/researchProject/integration/originals/human4_6W.RData")

# Load 7-11 weeks Seurat object
load("/Users/mayongzhi/Desktop/researchProject/integration/originals/human7_11W.RData")
human7_11W <- seurat
rm(seurat)

# Load 12-20 weeks Seurat object
load("/Users/mayongzhi/Desktop/researchProject/integration/originals/human12_20W.RData")
human12_20W <- humanpancreas.combined.sct
rm(humanpancreas.combined.sct)
ls()

# Explore the Seurat object
# Check what assays are available
Assays(human12_20W)

# View the metadata
head(human12_20W@meta.data)

# Access the RNA counts matrix
GetAssayData(human12_20W, assay = "RNA", slot = "counts")[1:5, 1:5]

# Check UMAP embeddings
Embeddings(human12_20W, "umap")[1:5, ]

# View cluster identities
Idents(human12_20W)

# Summary of dimensionality reduction
human12_20W[["pca"]]

table(human12_20W$seurat_clusters)
table(human12_20W@meta.data[["integrated_snn_res.0.1"]])
DimPlot(human12_20W, group.by = "Gestation_Age")

# See the number of variable feature in SCT assay
VariableFeatures(human12_20W[["SCT"]])[1:10]  # Show first 10 HVGs
length(VariableFeatures(human12_20W[["SCT"]])) 

# The number of variable feature in integrated assay
VariableFeatures(human12_20W[["integrated"]])[1:10]
length(VariableFeatures(human12_20W[["integrated"]]))

# The number of genes and cells before and after SCT
dim(GetAssayData(human12_20W, assay = "RNA", slot = "counts"))  # All genes
dim(GetAssayData(human12_20W, assay = "SCT", slot = "data"))     # Genes in SCT

# Check if SCT and integrated cover the same Gestation Age
DefaultAssay(human12_20W) <- "SCT"
table(human12_20W$Gestation_Age)

DefaultAssay(human12_20W) <- "integrated"
table(human12_20W$Gestation_Age)
