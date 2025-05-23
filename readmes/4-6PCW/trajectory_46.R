# 5- Trajectory analysis

# load dependencies
library(tidyverse)
library(Seurat)
#BiocManager::install("kstreet13/slingshot")
library(slingshot)
#BiocManager::install("destiny")
library(destiny)
#BiocManager::install("tradeSeq")
library(tradeSeq)

obj <- readRDS("~/Desktop/researchProject/integration/outputs/QC_seurat/human4_6W.rds")

DimPlot(obj, reduction = "umap", label = T) + NoLegend()

# filtering out cells
obj.filt <- subset(obj, seurat_clusters %in% c(5, 12, 19, 25, 28))
DimPlot(obj.filt, reduction = "umap", label = T)

# Generate the marker dict from paper
markers <- list(
  epithelial = c("PDX1", "EPCAM"),
  dorsal_MP = c("NR2F1"),
  ventral_MP = c("TBX3"),
  PB_progenitors = c("ISL1", "SULT1E1", "HHEX", "NKX6-2"),
  EHBD = c("SPP1"),
  enterocyte = c("CDX2"),
  hepatoblast = c("ALB")
)

obj.filt <- AddModuleScore(obj,
                           features = marker_dict,
                           name = "Score")

obj.filt@meta.data <- obj.filt@meta.data %>%
  dplyr::rename(Epithelial_Score = Score1,
                Mesenchymal_Score = Score2,
                Endothelial_Score = Score3,
                Neural_Score = Score4,
                Immune_Score = Score5,
                Erythroid_Score = Score6)

FeaturePlot(obj.filt, c("Epithelial_Score", "Mesenchymal_Score", "Endothelial_Score",
                        "Neural_Score", "Immune_Score", "Erythroid_Score"),
            label = TRUE, reduction = "umap",
            min.cutoff = 0)
#ggsave("~/Desktop/researchProject/integration/outputs/QC/46_celltypes_score.png")

# Annotate manually
new_names <- c("Mesenchymal","Mesenchymal","Mesenchymal","Epithelial","Mesenchymal","Mesenchymal","Mesenchymal","Mesenchymal","Mesenchymal",
  "Mesenchymal","Mesenchymal","Mesenchymal", "Mesenchymal",
  "Mesenchymal","Epithelial", "Endothelial","Mesenchymal",
  "Mesenchymal","Epithelial","Mesenchymal","Mesenchymal",
  "Epithelial","Mesenchymal","Mesenchymal","Neural",
  "Immune","Epithelial", "Epithelial", "Erythroid"
)

# Match new names to the seurat clusters
names(new_names) <- levels(obj.filt)
obj.filt <- RenameIdents(object = obj.filt, new_names)

DimPlot(obj.filt, label = TRUE) + NoLegend()
#ggsave("~/Desktop/researchProject/integration/outputs/QC/46_celltypeUMAP.png")

# Dot Plot for marker genes
marker_list <- c("ALB","AFP",
                 "SPP1","CDX2", "RFX6", 
                 "PDX1", "SOX9", "NKX6-1","EPCAM", 
                 "COL3A1","PECAM1", "PTPRC","ASCL1", "HBA1")
DotPlot(obj.filt, features = marker_list) + RotatedAxis()
#ggsave("~/Desktop/researchProject/integration/outputs/QC/46_markers.png")

DoHeatmap(obj.filt, features = marker_list)
colnames(obj.filt[[]])




## create a diffusion map
data <- GetAssayData(HB.seurat.filt, assay = "SCT", layer = "data")

data <- as.matrix(t(data))

data <- data[, VariableFeatures(HB.seurat.filt)]
dm <- DiffusionMap(data, n_pcs = 3)
HB.seurat.filt[["diffmap"]] <- CreateDimReducObject(embeddings = dm@eigenvectors,
                     key = "dc_",
                     assay = "SCT")

DimPlot(HB.seurat.filt, reduction = "diffmap", label = T)

# run slingshot
sds <- slingshot(Embeddings(HB.seurat.filt, "diffmap"),
          clusterLabels = HB.seurat.filt$seurat_clusters,
          start.clus = "4")

slingLineages(sds)

slingPseudotime(sds)

HB.seurat.filt[["pseudotime"]] <- slingPseudotime(sds)[,1]
FeaturePlot(HB.seurat.filt, "pseudotime",
            reduction = "umap", label = T)


FeatureScatter(HB.seurat.filt, "pseudotime", "STMN2") + 
  geom_smooth()

## tradeSeq for temporal gene comparison

pseudotime <- HB.seurat.filt$pseudotime[!is.na(HB.seurat.filt$pseudotime)]
weights <- slingCurveWeights(sds)[names(pseudotime),1]
counts <- GetAssayData(HB.seurat.filt, assay = "SCT", layer = "counts")
counts <- counts[,names(pseudotime)]

sce <- fitGAM(counts, pseudotime = pseudotime, cellWeights = weights)
sce
startRes <- startVsEndTest(sce)

# the top 10 temporally expressing genes
startRes %>% 
  filter(pvalue < 0.05) %>% 
  filter(waldStat > 100) %>% 
  slice_max(logFClineage1, n = 10)

startRes %>% 
  filter(pvalue < 0.05) %>% 
  filter(waldStat > 100) %>% 
  slice_min(logFClineage1, n = 10)

FeatureScatter(HB.seurat.filt, "pseudotime", "RBP1") +
  geom_smooth()

FeatureScatter(HB.seurat.filt, "pseudotime", "TTYH1") +
  geom_smooth()

saveRDS(HB.seurat, "~/Desktop/MSc_Applied_Bioinformatics/scRNAseq_workshop/tutorial/outputs/human_brain.rds")
saveRDS(HB.seurat.filt,"~/Desktop/MSc_Applied_Bioinformatics/scRNAseq_workshop/tutorial/outputs/human_brain_filt.rds" )
