# 4- Differential gene expression

# load packages
library(tidyverse)
library(Seurat)
library(patchwork)
options(future.globals.maxSize= 8 * 1024^3)
mem.maxVSize(vsize = 120000)

# Load Seurat object
obj <- readRDS("~/Desktop/researchProject/integration/outputs/QC_seurat/4_11PCW/human4_11W.rds")

# FindAllMarkers() was performed on HPC
#diff.markers <- FindAllMarkers(obj, test.use = "wilcox", only.pos = TRUE)
diff.markers <- read.csv("~/Desktop/researchProject/integration/metadata/markers_4_11W.csv")

FeaturePlot(obj, "PDX1", label = T)
FeaturePlot(obj, "SOX9")

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

DotPlot(obj, unique(unlist(markers))) + 
  theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
  axis.text.y = element_text(size = 8),
  plot.title = element_text(size = 10, face = "bold"))

# Dotplot for pancreatic marker genes
DotPlot(obj, c("PDX1", "NKX6-1", "SOX9")) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 10, face = "bold"))


FeaturePlot(obj, c("PDX1", "NKX6-1", "SOX9"), slot = "data", ncol = 3, 
            label = TRUE, label.size = 2)

# Feature plot for epithelial markers EPCAM and PDX1 
FeaturePlot(obj, c("EPCAM", "PDX1"), label = TRUE, label.size = 3)

# attach the diff.markers df to the Seurat object
Misc(obj, slot = "diff.markers") <- diff.markers
saveRDS(obj, "~/Desktop/researchProject/integration/outputs/QC_seurat/4_11PCW/human4_11W.rds")
