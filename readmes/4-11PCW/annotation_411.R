# Subcluster and annotate cell type

# load dependencies
library(tidyverse)
library(Seurat)
library(slingshot)
library(destiny)
library(tradeSeq)
library(harmony)
library(clustree)
library(EnhancedVolcano)
library(patchwork)
options(future.globals.maxSize= 8 * 1024^3)
mem.maxVSize(vsize = 120000)

obj <- readRDS("~/Desktop/researchProject/integration/outputs/QC_seurat/4_11PCW/human4_11W.rds")
DimPlot(obj, reduction = "umap", label = T) + NoLegend()

# Retain epithelial cells
obj.filt <- subset(obj, seurat_clusters %in% c(0, 5, 9, 10, 11, 13, 26, 29, 33))
DimPlot(obj.filt, reduction = "umap", label = T)

# PCA
obj.filt <- RunPCA(obj.filt, assay = "SCT", features = VariableFeatures(obj.filt))

# Harmony batch correction
obj.filt <- RunHarmony(object = obj.filt, group.by.vars = "orig.ident")

# ElbowPlot
ElbowPlot(obj.filt, ndims = 50)

# UMAP and clustering with Harmony embeddings
obj.filt <- RunUMAP(obj.filt, reduction = "harmony", dims = 1:50)
obj.filt <- FindNeighbors(obj.filt, reduction = "harmony", dims = 1:50)

# Make clustree
obj.filt@meta.data <- obj.filt@meta.data[, !grepl("^SCT_snn_res\\.", colnames(obj.filt@meta.data))]

res.range <- seq(0.8, 1.8,0.1)
res.range

for (res in res.range){
  obj.filt <- FindClusters(obj.filt, resolution = res)
}

# Visualise clusters using clustree
clustree(obj.filt, prefix = "SCT_snn_res.")

# Find cluster
obj.filt <- FindClusters(obj.filt, resolution = 1.2)

# Save the object
saveRDS(obj.filt, "~/Desktop/researchProject/integration/outputs/QC_seurat/4_11PCW/epith_411.rds")

# Find markers
# FindAllMarkers() was performed on HPC
diff.markers <- read.csv("~/Desktop/researchProject/integration/metadata/markers_epith_411/FindAllMarkers_results.csv")
Misc(obj.filt, slot = "diff.markers") <- diff.markers

top <- diff.markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 1)

# The marker dict from paper (We don't use it for annotation as there are too many markers)
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

# The markers Ma used for annotation
ma_markers <- c("PDX1", "NKX6-1", "PTF1A", "GATA4", "FOXA2", "TBX3", "NR2F1", 
                "CPA2", "RBPJL", "CTRB2", "CLPS", "CTRB1", "AMY2B", "HNF1B",
                "HES4", "DCDC2", "CFTR", "ASCL2", "NEUROG3", "NEUROD1", "INS",
                "GCG", "PPY", "SST", "GHRL")

early_epith_markers <- c("PDX1", "PTF1A","NKX6-1","NR2F1", "SOX6", "TBX3", "ISL1",
                         "NKX6-2", 'HHEX', "SULT1E1", "SPP1", "LGALS3", "CDX2", "RFX6", "ALB","APOA2")

combined_markers <- unique(c(ma_markers, early_epith_markers))

PB_markers <- c("PDX1", "HHEX", "ISL1", "SULT1E1", "NKX6-2")

celltype_markers <- list(
  dorsal_MP = c("PDX1", "NR2F1"),
  ventral_MP = c("PDX1", "TBX3"),
  PB_progenitors = c("HHEX","ISL1", "NKX6-2", "SULT1E1"),
  EHBD = c("SPP1", "SULT1E1"),
  enterocyte = c("CDX2"),
  hepatoblast = c("ALB")
)
obj.filt <- AddModuleScore(obj.filt,
                            features = celltype_markers,
                            name = "Score")

obj.filt@meta.data <- obj.filt@meta.data %>%  
  dplyr::rename(dorsalMP_Score = Score1,
                ventralMP_Score = Score2,
                PB_Score = Score3,
                EHBD_Score = Score4,
                Enterocyte_Score = Score5,
                Hepatoblast_Score = Score6,
                )

score_list <- c("dorsalMP_Score", "ventralMP_Score", "PB_Score", "EHBD_Score", "Enterocyte_Score", "Hepatoblast_Score")

# Generate individual violin plots with median dot and reference line
violin_plots <- lapply(score_list, function(score) {
  VlnPlot(obj.filt, features = score, pt.size = 0) +
    stat_summary(fun = median, geom = "point", size = 2, colour = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    NoLegend() +
    ggtitle(score)
})

# Combine plots in 2 rows
combined_plot <- wrap_plots(violin_plots, ncol = 2)
combined_plot

DotPlot(obj.filt, combined_markers, group.by = "seurat_clusters", dot.scale = 5) + 
  scale_color_gradient (low = "darkslateblue", high = "brown2")+ 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 10, face = "bold"))

DotPlot(obj.filt, ma_markers, dot.scale = 4) + 
  scale_color_gradient (low = "darkslateblue", high = "brown2")+ 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 10, face = "bold"))

DotPlot(obj.filt, early_epith_markers, dot.scale = 4) + 
  scale_color_gradient (low = "darkslateblue", high = "brown2")+ 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 10, face = "bold"))

DotPlot(obj.filt, c("SST", "GHRL"))

# Annotate manually
Idents(obj.filt) <- "seurat_clusters"

# Subset to keep only cells whose cluster ID is NOT “11”
obj.filt <- subset(obj.filt, idents = setdiff(levels(Idents(obj.filt)), "11"))

DimPlot(obj.filt, label = T)

new_names <- c("Early_tip", "Early_trunk", "Trunk", "EHBD", "ventral_MP",
               "Acinar", "Early_tip", "Tip", "Early_tip", "Early_tip", 
               "Early_tip", "PB", "Early_trunk", "EP", "Early_tip",
               "Early_tip", "Duct", "Tip", "Hepatoblast", "Early_tip", "Beta", 
               "Early_trunk", "dorsal_MP", "Alpha"
)

# Match new names to the seurat clusters
names(new_names) <- levels(obj.filt)
obj.filt <- RenameIdents(object = obj.filt, new_names)
DimPlot(obj.filt, label = TRUE) + NoLegend()


# Dot Plot for marker genes
Idents(obj.filt)<- factor(Idents(obj.filt), levels = c("PB", "EHBD", "ventral_MP", "dorsal_MP",
                            "Early_tip", "Tip", "Acinar",
                            "Early_trunk", "Trunk", "Duct",
                            "EP", "Alpha", "Beta", "Hepatoblast"))



DotPlot(obj.filt, c("SOX9", combined_markers), dot.scale = 5) + 
  scale_color_gradient (low = "darkslateblue", high = "brown2")+ 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 10, face = "bold"))

# Save the object
saveRDS(obj.filt, "~/Desktop/researchProject/integration/outputs/QC_seurat/4_11PCW/epith_411.rds")

