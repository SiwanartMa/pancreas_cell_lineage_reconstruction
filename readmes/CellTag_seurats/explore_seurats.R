library(Seurat)
library(glue)
library(patchwork)
library(BPCells)
library(scDblFinder)
library(tidyverse)
library(dplyr)
#install.packages('BPCells', repos = c('https://bnprks.r-universe.dev', 'https://cloud.r-project.org'))


# import the celltagged Seurat objects
for (i in c(2, 6, 8, 10)) {
  path <- glue("/Users/mayongzhi/Desktop/researchProject/integration/originals/ScaleRNA/CellTag-D{i}.ScaleRNA_SeuratObject.rds")
  obj <- readRDS(path)
  assign(glue("d{i}_tag"), obj)
}

ctrl_tag <- readRDS("/Users/mayongzhi/Desktop/researchProject/integration/originals/ScaleRNA/CellTag-iPSCs-Frozen.ScaleRNA_SeuratObject.rds")

# Merge all objects
obj <- merge(x = ctrl_tag, y = list(d2_tag, d6_tag, d8_tag, d10_tag))

# Assume `obj` is your Seurat object
metadata <- as.data.frame(obj@meta.data)
class(metadata)

# Create mapping from each RT barcode (sample) to its ligation barcodes
rt_ligation_map <- metadata %>%
  select(RT, Ligation) %>%
  group_by(RT) %>%
  dplyr::summarise(
    Ligation_barcodes = paste(unique(Ligation), collapse = ", "),
    num_cells = n(),
    .groups = "drop"
  )

# View the mapping
print(rt_ligation_map)

# Create a mapping from i5 wells to their RT barcodes
i5_mapping <- metadata %>%
  select(RT, Ligation, i5) %>%
  group_by(i5) %>%
  summarise(
    RT_barcodes = paste(unique(RT), collapse = ", "),
    num_cells = n(),
    .groups = "drop"
  )

# View the mapping
print(i5_mapping)


