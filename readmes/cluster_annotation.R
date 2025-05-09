# Install packages
library(Seurat)
library(tidyverse)
library(patchwork)
#install.packages("scCustomize")
library(scCustomize)
#install.packages("qs")
library(qs)


# Load the 3 Seurat objects
# Increase the memory
mem.maxVSize(vsize = 120000)

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

human4_6W
table(human4_6W@meta.data$orig.ident)
markers <- FindAllMarkers(human4_6W, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
str(markers)
all_markers <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  arrange(cluster)

top_5 <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 5, data_frame = TRUE,rank_by = "avg_log2FC")
head(top_5, 10)

# Read in pre-computed data.frame results
presto_markers <- qread("assets/presto_markers.qs")

# Extract top markers
top_5_presto <- Extract_Top_Markers(marker_dataframe = presto_markers, num_genes = 5, group_by = "cluster",
                                    rank_by = "auc", gene_column = "gene")
head(top_5_presto, 10)
Iterate_FeaturePlot_scCustom(seurat_object = human4_6W, gene_list = top_5, single_pdf = F)
print(top_markers)
#write.csv(markers, file = "/Users/mayongzhi/Desktop/researchProject/integration/outputs/4_6W_cluster_markers.csv", row.names = FALSE)
cluster_ids <- levels(as.factor(markers$cluster))  
new_ids <- c(
  "Endothelial",    # cluster 0: e.g., PECAM1, SPARCL1
  "Mesenchymal",    # cluster 1: e.g., COL3A1, COL6A3
  "Endocrine",      # cluster 2: CHGA, NEUROD1
  "Immune",         # cluster 3: PTPRC
  "Erythroid",      # cluster 4: HBA1
  "Neural",         # cluster 5: ASCL1
  "Non-endocrine",  # cluster 6: CPA2, PDX1
  "Other"           # cluster 7: unclassified or unknown
)

# Map to Seurat object
names(new_ids) <- cluster_ids
human4_6W <- RenameIdents(human4_6W, new_ids)

# Save annotation to metadata
human4_6W$celltype <- Idents(human4_6W)

# Plot annotated UMAP
DimPlot(human4_6W, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend()
genes <- c("CPA2", "PDX1", "CHGA","NEUROD1", "EPCAM", "COL3A1", "PECAM1", "PTPRC", "ASCL1", "HBA1")
DotPlot(object = human4_6W, features=genes, group.by = "days")
FeaturePlot(human4_6W,
            features = genes,
            label.size=2)
