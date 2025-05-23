# load packages
library(tidyverse)
library(Seurat)
library(patchwork)
library(harmony)
library(clustree)
options(future.globals.maxSize= 8 * 1024^3)
mem.maxVSize(vsize = 80000)

# Load Seurat object
obj <- readRDS("~/Desktop/researchProject/integration/outputs/QC_seurat/human4_6W.rds")

# Subset for PCW 4–5 epithelial cells
epith_4_5 <- subset(obj, subset = days %in% c("W4", "W5") & seurat_clusters %in% c(5,12,19, 25, 28))

# Normalize and find variable features
epith_4_5 <- NormalizeData(epith_4_5)
epith_4_5 <- FindVariableFeatures(epith_4_5, selection.method = "vst", nfeatures = 2000)

# Scale data and regress out unwanted variables
epith_4_5 <- ScaleData(epith_4_5, vars.to.regress = c("S.Score", "G2M.Score", "nCount_RNA"))

# PCA
epith_4_5 <- RunPCA(epith_4_5, features = VariableFeatures(epith_4_5))

# Harmony batch correction
epith_4_5 <- RunHarmony(object = epith_4_5, 
                        group.by.vars = "orig.ident")

# UMAP and clustering with Harmony embeddings
epith_4_5 <- RunUMAP(epith_4_5, reduction = "harmony", dims = 1:30)
epith_4_5 <- FindNeighbors(epith_4_5, reduction = "harmony", dims = 1:30)

# Make clustree
res.range <- seq(0.1, 1,0.1)
res.range

for (res in res.range){
  epith_4_5 <- FindClusters(epith_4_5, resolution = res)
}

# Visualise clusters using clustree
clustree <- clustree(epith_4_5, prefix = "SCT_snn_res.")
clustree
#ggsave("~/Desktop/researchProject/integration/outputs/QC_plot/4_6PCW_reQC/epith_4_5_clustree.png")

epith_4_5 <- FindClusters(epith_4_5, resolution = 0.7)

# The markers from paper
markers <- list(
  dorsal_MP = c("PDX1", "NR2F1"),
  ventral_MP = c("PDX1","TBX3"),
  PB_progenitors = c("ISL1", "SULT1E1"),
  EHBD = c("SPP1", "SULT1E1"),
  enterocyte = c("CDX2"),
  hepatoblast = c("ALB")
)

epith_4_5 <- AddModuleScore(epith_4_5,
                            features = markers,
                            name = "Score")

epith_4_5[[]]

epith_4_5@meta.data <- epith_4_5@meta.data %>%  
  dplyr::rename(dorsalMP_Score = Score1,
                ventralMP_Score = Score2,
                PB_Score = Score3,
                EHBD_Score = Score4,
                Enterocyte_Score = Score5,
                Hepatoblast_Score = Score6)

FeaturePlot(epith_4_5, c("PB_Score", "dorsalMP_Score", "ventralMP_Score", "EHBD_Score", "Enterocyte_Score", "Hepatoblast_Score"), label = TRUE, 
            min.cutoff = 0, reduction = "umap")
#ggsave("~/Desktop/researchProject/integration/outputs/QC_plot/4_6PCW_reQC/epith_4_5_feature.png")

# List of all scores to plot
score_list <- c("dorsalMP_Score", "ventralMP_Score", "PB_Score", "EHBD_Score", "Enterocyte_Score", "Hepatoblast_Score")

# Generate individual violin plots with median dot and reference line
violin_plots <- lapply(score_list, function(score) {
  VlnPlot(epith_4_5, features = score, pt.size = 0) +
    stat_summary(fun = median, geom = "point", size = 2, colour = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    NoLegend() +
    ggtitle(score)
})

# Combine plots in 2 rows
combined_plot <- wrap_plots(violin_plots, ncol = 2)
#ggsave("~/Desktop/researchProject/integration/outputs/QC_plot/4_6PCW_reQC/epith_4_5_viol.png")

# Display the plot
print(combined_plot)


# Annotate clusters
new.cluster.ids <- c("Enterocyte", "dorsalMP", "ventralMP", "Enterocyte", "PB + EHBD",
                     "EHBD", "EHBD", "ventralMP", "Enterocyte", "Enterocyte")

names(new.cluster.ids) <- levels(epith_4_5)
epith_4_5 <- RenameIdents(epith_4_5, new.cluster.ids)
DimPlot(epith_4_5, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#ggsave("~/Desktop/researchProject/integration/outputs/QC_plot/4_6PCW_reQC/epith_4_5_celltype.png")

# Pseudotime
data <- as.matrix(t(epith_4_5@assays$SCT@data))
data <- data[, VariableFeatures(epith_4_5)]
dm <- DiffusionMap(data, n_pcs = 3)

epith_4_5[["diffmap"]] <- CreateDimReducObject(embeddings = dm@eigenvectors[,1:3], key = "dc_", assay = "SCT")
DimPlot(epith_4_5, reduction = "diffmap", label = TRUE) 

sds <- slingshot(Embeddings(epith_4_5, "diffmap"), 
                 clusterLabels = epith_4_5$seurat_clusters,
                 start.clus= '4')

slingLineages(sds)
slingPseudotime(sds)[,1]

epith_4_5$pseudotime <- slingPseudotime(sds)[,1]

FeaturePlot(epith_4_5, "pseudotime", reduction = "umap", label = T,) +
  scale_color_gradient(low = "lightgreen", high = "orange")
#ggsave("~/Desktop/researchProject/integration/outputs/QC_plot/4_6PCW_reQC/epith_4_5_pseudotime.png")

FeatureScatter(epith_4_5, "pseudotime", "ISL1") + 
  geom_smooth() +
  ggtitle("ISL1 Expression vs. Pseudotime")

FeatureScatter(epith_4_5, "pseudotime", "SULT1E1") + 
  geom_smooth()+
  ggtitle("SULT1E1 Expression vs. Pseudotime")

# Save processed object
saveRDS(epith_4_5, file = "~/Desktop/researchProject/integration/outputs/QC_seurat/epithelial_PWC4_5_harmony.rds")

