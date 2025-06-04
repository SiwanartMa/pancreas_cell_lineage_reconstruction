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
library(harmony)
library(clustree)

obj <- readRDS("~/Desktop/researchProject/integration/outputs/QC_seurat/4_11PCW/human4_11W.rds")
DimPlot(obj, reduction = "umap", label = T) + NoLegend()

# Retain epithelial cells
obj.filt <- subset(obj, seurat_clusters %in% c(0, 5, 9, 10, 11, 13, 26, 29))
DimPlot(obj.filt, reduction = "umap", label = T)
DimPlot(obj.filt, reduction = "umap", group.by = "days", label = T)

# Normalize and find variable features
obj.filt <- NormalizeData(obj.filt)
obj.filt <- FindVariableFeatures(obj.filt, selection.method = "vst", nfeatures = 3000)

# Scale data and regress out unwanted variables
obj.filt <- ScaleData(obj.filt, vars.to.regress = c("S.Score", "G2M.Score", "nCount_RNA"))

# PCA
obj.filt <- RunPCA(obj.filt, features = VariableFeatures(obj.filt))

# Harmony batch correction
obj.filt <- RunHarmony(object = obj.filt, 
                        group.by.vars = "orig.ident")

# ElbowPlot
ElbowPlot(obj.filt, ndims = 50)

# UMAP and clustering with Harmony embeddings
obj.filt <- RunUMAP(obj.filt, reduction = "harmony", dims = 1:50)
obj.filt <- FindNeighbors(obj.filt, reduction = "harmony", dims = 1:50)

# Make clustree
obj.filt@meta.data <- obj.filt@meta.data[, !grepl("^SCT_snn_res\\.", colnames(obj.filt@meta.data))]

res.range <- seq(1.2, 2.1,0.1)
res.range

for (res in res.range){
  obj.filt <- FindClusters(obj.filt, resolution = res)
}

# Visualise clusters using clustree
clustree(obj.filt, prefix = "SCT_snn_res.")

# Find cluster
obj.filt <- FindClusters(obj.filt, resolution = 1.7)

# Save the object
saveRDS(obj.filt, "~/Desktop/researchProject/integration/outputs/QC_seurat/4_11PCW/epith_411.rds")

# Find markers
# FindAllMarkers() was performed on HPC
diff.markers <- read.csv("~/Desktop/researchProject/integration/metadata/markers_epith_411.csv")
Misc(obj.filt, slot = "diff.markers") <- diff.markers

top <- diff.markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 1)
  
# Dimplot
DimPlot(obj.filt, label = T) + NoLegend()

# Generate the marker dict from paper
markers<- list(
  Epsilon = c("ARX", "ETV1", "INSM1", "ISL1", "NEUROD1", "NKX2-2", "RFX6", "ST18"),
  EP = c("NEUROD1", "INSM1", "NKX2-2", "HES6", "ISL1", "MAFB", "PAX4", "SOX4", 
         "PAX6", "FEV", "RFX6","MLXIPL", "DACH1", "NEUROG3", "TOX"),
  Beta = c("MAFB", "PLAGL1", "PAX6", "ISL1", "NEUROD1", "INSM1", "MEIS2", 
           "NKX2-2", "MAFA", "MLXIPL"),
  Delta = c("ISL1", "NEUROD1", "MAFB", "INSM1", "HHEX", "NKX2-2", "FOXO1", "POU2F2"),
  Alpha = c("ISL1", "MAFB", "PAX6", "NEUROD1", "IRX2", "ARX", "INSM1", "FEV", 
            "MLXIPL", "NKX2-2"),
  Duct = c("ID3", "ID1", "HES4", "RBPJ", "HEY1", "HIF3A", "FOSB", "ID4", "NR4A1", 
           "EHF", "ELF3", "FOS", "KLF6", "NFIB"),
  Trunk = c("ASCL2", "ID4", "ELF3", "SOX4", "HES1", "HES4", "ATF3", "NFIB", "EHF", 
            "GLIS3", "ID1"),
  Dorsal_MP = c("NR2F1", "ID3", "HES1", "ZFP36L2", "ID1", "STAT1", "ONECUT2", 
                "ZFP36L1", "SOX11", "SALL4", "FOXA1", "ID2", "MEIS1", "PDX1", 
                "HOXB2", "FOXA2", "SOX9", "ONECUT1", "HMGA1", "GATA4"),
  Ventral_MP = c("ID3", "ZFP36L2", "ID1", "SOX11", "ID2", "STAT1", "PTF1A", "ZFP36L1", 
                 "TBX3", "HES1", "PDX1", "ONECUT3", "ONECUT1", "SOX9", "GATA4", 
                 "NKX6-1", "HMGA1", "FOXA2", "SALL4", "LIN28A"), 
  Early_trunk = c("HMGB2", "HES4", "HMGB1", "HMGA1", "HMGB3"),
  Early_tip = c("HMGA1", "NFE2L3", "HMGB3", "YBX1"),
  Tip = c("FOSB", "JUNB", "FOS", "NR4A1", "ATF3", "JUN", "MAFF", "RBPJ", 
          "EGR1", "YBX3", "RBPJ", "KLF6", "JUND", "CEBPB", "SOX6", "XBP1", 
          "EGR3", "STAT3", "EGR2", "ZNF385B", "EPAS1", "PEG3", "KLF10"), 
  Acinar = c("ATF3", "FOSB", "JUNB", "KLF6", "RBPJ", "XBP1", "MAFF", "FOS", 
             "NR4A1", "MYC", "JUND", "EGR1", "CSRNP1")
)

obj.filt <- AddModuleScore(obj.filt,
                           features = markers,
                           name = "Score")


obj.filt@meta.data <- obj.filt@meta.data %>%
  dplyr::rename(Epsilon_Score = Score1,
                EP_Score = Score2,
                Beta_Score = Score3,
                Delta_Score = Score4,
                Alpha_Score = Score5,
                Duct_Score = Score6, 
                Trunk_Score = Score7,
                DorsalMP_Score = Score8,
                VentralMP_Score = Score9, 
                EarlyTrunk_Score = Score10,
                EarlyTip_Score = Score11,
                Tip_Score = Score12,
                Acinar_Score = Score13)

celltype.list <- c("Epsilon_Score","EP_Score","Beta_Score", "Delta_Score", "Alpha_Score",
                   "Duct_Score", "Trunk_Score", "DorsalMP_Score", "VentralMP_Score",
                   "EarlyTrunk_Score", "EarlyTip_Score", "Tip_Score", "Acinar_Score")


Magdy_markers <- c("ISL1", "CHGA", "CHGB", "PAX6", 
                   "NGN3", "FEV", 
                   "GCG",
                   "INS", 
                   "SPP1", "SOX9", "PROX1", "ONECUT1",
                   "MKI67", 
                   "CPA1", "CPA2", "MECOM", "GATA4")

Ma_markers <- c("PDX1", "NKX6-1", "PTF1A", "GATA4", "FOXA2", "TBX3", "NR2F1", 
                "CPA2", "RBPJL", "CTRB2", "CLPS", "CTRB1", "AMY2B", "HNF1B",
                "HES4", "DCDC2", "CFTR", "ASCL2", "NEUROG3", "NEUROD1", "INS",
                "GCG", "PPY", "SST", "GHRL")

Combined_markers <- c(unique(unlist(markers)), unique(Magdy_markers))

DotPlot(obj.filt, Ma_markers, dot.scale = 4) + 
  scale_color_gradient (low = "darkslateblue", high = "brown2")+ 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 10, face = "bold"))



# Annotate manually
obj.filt <- subset(obj.filt, seurat_clusters != 19)
DimPlot(obj.filt, label = T)

new_names <- c("Early_trunk", "Early_tip", "Trunk", "Early_tip", "Early_tip",
               "Early_trunk", "Tip", "Dorsal_MP", "Early_tip", "Acinar", "Early_trunk",
               "Early_trunk", "Trunk", "Early_tip", "Acinar", "Tip",
               "Ventral_MP", "Acinar", "EP", "Beta", "Tip", "Early_tip", 
               "Alpha", "Delta", "Epsilon"
)

# Match new names to the seurat clusters
names(new_names) <- levels(obj.filt)
obj.filt <- RenameIdents(object = obj.filt, new_names)

DimPlot(obj.filt, label = TRUE) + NoLegend()


# Dot Plot for marker genes
DotPlot(obj.filt, Ma_markers, dot.scale = 5) + 
  scale_color_gradient (low = "darkslateblue", high = "brown2")+ 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 10, face = "bold"))


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
