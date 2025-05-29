# Seurat clustering
library(tidyverse)
library(Seurat)
library(clustree)

# Load Seurat object
obj <- readRDS("~/Desktop/researchProject/integration/outputs/ScaleRNA/QC_seurat/combined.rds")

# plot a scatter of nCount_RNA vs TUBB
FeatureScatter(obj, "nCount_RNA", "TUBB", slot = "counts")
# It is positive correlate with the sequencing depth

# run SCTransform workflow
obj <- SCTransform(obj, method = "glmGamPoi")

# PCA
obj <- RunPCA(obj)
# extract out the best 50 dimensions and save it in a new slot

# visualise dim reduction
DimPlot(obj)
DimPlot(obj, group.by = "sample")

# Colourcode cells by gene expression
FeaturePlot(obj, "TUBB")
FeaturePlot(obj, "percent_mito") 
FeaturePlot(obj, "G2M.Score")

# rerun SCTransform, with regression
obj <- SCTransform(obj,
                   method = "glmGamPoi",
                   vars.to.regress = c("nCount_RNA", "S.Score", "G2M.Score",
                                       "percent_mito", "percent_ribo"))

# rerun PCA
obj <- RunPCA(obj)
DimPlot(obj, group.by = "sample") # Cells in different phases overlap
FeaturePlot(obj, "percent_mito") # The bias is better 

## Non-linear dimensional reduction
obj <- RunUMAP(obj, dims = 1:50)
obj <- RunTSNE(obj, dims = 1:50)

# Plot UMAP 
DimPlot(obj, group.by = "sample", reduction = "umap")

# Add 'group' column based on sample string content
obj$group <- ifelse(grepl("CellTag", obj$sample), "CellTag", "CTRL")

# Add 'day' column based on sample string content
obj$day <- case_when(
  grepl("D2", obj$sample, ignore.case = TRUE) ~ "D2",
  grepl("D6", obj$sample, ignore.case = TRUE) ~ "D6",
  grepl("D8", obj$sample, ignore.case = TRUE) ~ "D8",
  grepl("D10", obj$sample, ignore.case = TRUE) ~ "D10",
  grepl("iPSCs-Frozen", obj$sample, ignore.case = TRUE) ~ "iPSC",
  TRUE ~ "Other"
)

DimPlot(obj, group.by = "day", reduction = "umap", split.by = "group")

### Clustering
obj <- FindNeighbors(obj, dims = 1:50)

res.range <- seq(0.1, 1.5,0.1)
res.range

for (res in res.range){
  obj <- FindClusters(obj, resolution = res)
}

head(obj[[]])

# Visualise clusters using clustree
clustree(obj, prefix = "SCT_snn_res.") 

# cluster using res.1.6
obj <- FindClusters(obj, resolution=1.4)

DimPlot(obj, group.by = "group", reduction = "umap", label = F)

saveRDS(obj, "~/Desktop/researchProject/integration/outputs/ScaleRNA/QC_seurat/combined.rds")
