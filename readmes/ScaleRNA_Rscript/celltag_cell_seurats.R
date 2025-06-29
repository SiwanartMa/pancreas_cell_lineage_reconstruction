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

# Set working directory
setwd("~/Desktop/researchProject/integration/originals/ScaleRNA/")

# Import celltag_count.csv
celltag_counts <- read_csv("~/Desktop/researchProject/integration/metadata/celltag_counts.csv")
celltag_counts$Cell.BC

obj <- obj.filt

# Extract just the barcode part from obj's column names
clean_barcodes <- gsub("_CellTag-[^_]+\\.ScaleRNA$", "", colnames(obj))

# Create logical vector: TRUE if stripped barcode matches celltag_counts$Cell.BC
obj$express_celltag <- clean_barcodes %in% celltag_counts$Cell.BC

table(obj$express_celltag)

obj_celltag <- subset(obj, express_celltag == TRUE)
rownames(obj_celltag[[]])


DimPlot(
  obj,
  group.by = "seurat_clusters",
  label = TRUE,
  cells.highlight = colnames(obj)[obj$express_celltag],
  sizes.highlight = 2,
  cols.highlight = "red",
  cols = "lightgray"
)


