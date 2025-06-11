library(Seurat)
library(tidyverse)
library(patchwork)
library(devtools)
#install_github("UPSUTER/GEMLI", subdir="GEMLI_package_v0")
library(GEMLI)
library(igraph)
library(UpSetR)
# Increase the memory
mem.maxVSize(vsize = 200000)

# Load 4-11 weeks Seurat object
obj <- readRDS("/Users/mayongzhi/Desktop/researchProject/integration/outputs/QC_seurat/4_11PCW/epith_411.rds")
GEMLI_items <- readRDS("/Users/mayongzhi/Desktop/researchProject/integration/outputs/gemli/epith_411_gemli.rds")
head(GetAssayData(obj))
GEMLI_items$gene_expression[1:5]

#gene_expression_matrix <- as.matrix(GetAssayData(obj, assay = "SCT", layer = "data"))
#gene_expression_matrix[1:5, 1:5]
#dim(gene_expression_matrix)

#GEMLI_items <- list(
#  gene_expression = gene_expression_matrix
#)

# Predict lineages
#GEMLI_items = predict_lineages(GEMLI_items)
#GEMLI_items[['prediction']][1:5,15:19]

# Make celltype dataframe
celltype.df <- data.frame(
  cell.ID = colnames(obj), 
  cell.type = obj@active.ident,
  row.names = NULL
)

GEMLI_items[['cell_type']] <- celltype.df

# Assign cell type color
cell.type <- unique(GEMLI_items[['cell_type']]$cell.type)
color <- c("#5386BD", "skyblue1", "darkgreen", "gold", "red", "purple", "brown",
           "orange", "coral2", "deepskyblue4", "cyan4", "aquamarine", "azure3", "chocolate")
Cell_type_color <- data.frame(cell.type, color)
GEMLI_items[['cell_type_color']] = Cell_type_color

# Visualize as network
visualize_as_network(GEMLI_items, cutoff=95, display_orphan=F, max_edge_width=1, ground_truth=F, include_labels=F, layout_style="fr", cell_type_colors = T)
visualize_as_network(GEMLI_items, cutoff=99, display_orphan=F, max_edge_width=1, ground_truth=F, include_labels=F, layout_style="kk", cell_type_colors = T)
visualize_as_network(GEMLI_items, cutoff=99, display_orphan=F, max_edge_width=1, ground_truth=F, include_labels=T, layout_style="grid", cell_type_colors = T)

# Plot cell type composition
GEMLI_items = prediction_to_lineage_information(GEMLI_items, cutoff=50)
cell_type_composition_plot(GEMLI_items, ground_truth=F, cell_type_colors=T, type=c("upsetR")) 

# Extrace cell fate lineages (Not success yet)
GEMLI_items
GEMLI_items<-extract_cell_fate_lineages(GEMLI_items, selection=c("EP", "Beta"), unique=FALSE, threshold = c(10,10))
GEMLI_items[['cell_fate_analysis']][1:10,] 
table(GEMLI_items[['cell_fate_analysis']]$cell.fate)

Lookup <- merge(GEMLI_items[['predicted_lineage_table']], GEMLI_items[['cell_type']], by = "cell.ID", all = TRUE)
Lookup

GEMLI_items$cell_fate_analysis
GEMLI_items<-cell_fate_DEG_calling(GEMLI_items, ident1="asym_EP", ident2="sym_EP", layer = "counts", min.pct=0.05, logfc.threshold=0.1)
GEMLI_items[['DEG']][1:10,]
