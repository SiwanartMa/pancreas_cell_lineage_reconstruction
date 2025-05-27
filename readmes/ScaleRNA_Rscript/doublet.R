library(Seurat)
library(glue)
library(patchwork)
library(BPCells)
library(scDblFinder)
library(tidyverse)
library(dplyr)


# import the celltagged Seurat objects
for (i in c(2, 6, 8, 10)) {
  path <- glue("/Users/mayongzhi/Desktop/researchProject/integration/originals/ScaleRNA/CellTag-D{i}.ScaleRNA_SeuratObject.rds")
  obj <- readRDS(path)
  assign(glue("d{i}_tag"), obj)
}

ctrl_tag <- readRDS("/Users/mayongzhi/Desktop/researchProject/integration/originals/ScaleRNA/CellTag-iPSCs-Frozen.ScaleRNA_SeuratObject.rds")

# Get count matrix
setwd("/Users/mayongzhi/Desktop/researchProject/integration/originals/ScaleRNA/")

# Input list of Seurat objects
obj.names <- c("d2_tag", "d6_tag", "d8_tag", "d10_tag", "ctrl_tag")

for (name in obj.names) {
  # Get the object by name
  obj <- get(name)
  
  # Convert counts to dense matrix
  raw_counts <- as.matrix(GetAssayData(obj, layer = "counts"))
  
  # Run scDblFinder
  sce <- scDblFinder(raw_counts)
  obj <- AddMetaData(obj,
                     metadata = sce$scDblFinder.score,
                     col.name = "scDblFinder.score")
  obj <- AddMetaData(obj,
                     metadata = sce$scDblFinder.class,
                     col.name = "scDblFinder.class")
  
  # Print classification table
  cat("Processed:", name, "\n")
  print(table(obj[["scDblFinder.class"]]))
  
  # Update object in environment
  assign(name, obj, envir = .GlobalEnv)
}

# Save objects
output_dir <- "~/Desktop/researchProject/integration/outputs/ScaleRNA/QC_seurat/"

for (name in obj.names){
  obj <- get(name)
  saveRDS(obj, file = file.path(output_dir, paste0(name, ".rds")))
}


