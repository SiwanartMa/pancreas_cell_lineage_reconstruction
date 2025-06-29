library(Seurat)
library(tidyverse)
library(patchwork)
library(devtools)
#install_github("UPSUTER/GEMLI", subdir="GEMLI_package_v0")
library(GEMLI)
library(igraph)
library(UpSetR)
library(ggrepel)
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

# subset GEMLI by cell types
# Extract the cell IDs of interest
tip_cells <- GEMLI_items$cell_type$cell.ID[GEMLI_items$cell_type$cell.type %in% c("Early_tip", "Tip", "Acinar")]


# Subset all relevant parts of the list using these cell IDs
GEMLI_tip <- list(
  gene_expression = GEMLI_items$gene_expression[, tip_cells, drop = FALSE],
  prediction = GEMLI_items$prediction[tip_cells, tip_cells, drop = FALSE],
  cell_type = GEMLI_items$cell_type[GEMLI_items$cell_type$cell.ID %in% tip_cells, ],
  cell_type_color = GEMLI_items$cell_type_color,
  predicted_lineages = GEMLI_items$predicted_lineages[tip_cells],
  predicted_lineage_table = GEMLI_items$predicted_lineage_table[GEMLI_items$predicted_lineage_table[, "cell.ID"] %in% tip_cells, ]
)

GEMLI_tip <- prediction_to_lineage_information(GEMLI_tip, cutoff=95)
cell_type_composition_plot(GEMLI_tip, ground_truth=F, cell_type_colors=T, type=c("upsetR")) 


# Plot cell type composition
GEMLI_items = prediction_to_lineage_information(GEMLI_items, cutoff=95)
cell_type_composition_plot(GEMLI_items, ground_truth=F, cell_type_colors=T, type=c("upsetR")) 


# Extrace cell fate lineages
GEMLI_tip<-extract_cell_fate_lineages(GEMLI_tip, selection=c("Tip", "Acinar"), unique=FALSE, threshold = c(10,10))
GEMLI_tip[['cell_fate_analysis']][1:10,] 
table(GEMLI_tip[['cell_fate_analysis']]$cell.type)

GEMLI_tip<-cell_fate_DEG_calling(GEMLI_tip, ident1="sym_Tip", ident2="asym_Tip", min.pct=0.05, logfc.threshold=0.1)
GEMLI_tip[['DEG']][1:10,]


# Perform cell fate DEG calling
ident1 <- "sym_Tip"
ident2 <- "asym_Tip"
min.pct <- 0.05
logfc.threshold <- 0.1

GEMLI_Seurat<-CreateSeuratObject(GEMLI_tip[['gene_expression']], project = "SeuratProject", assay = "RNA")
Metadata<-GEMLI_tip[['cell_fate_analysis']];Metadata$ident<-NA
Metadata$ident[Metadata$cell.fate %in% ident1]<-"ident1"
Metadata$ident[Metadata$cell.fate%in%ident2]<-"ident2"
Meta<-as.data.frame(Metadata[,c(5)]); rownames(Meta)<-Metadata$cell.ID; colnames(Meta)<-c("cell.fate")
GEMLI_Seurat<-AddMetaData(GEMLI_Seurat, Meta, col.name = NULL); DefaultAssay(object = GEMLI_Seurat) <- "RNA"; Idents(GEMLI_Seurat) <- GEMLI_Seurat$cell.fate
GEMLI_Seurat <- NormalizeData(GEMLI_Seurat)
DEG <- FindMarkers(object = GEMLI_Seurat, ident.1 = "ident1", ident.2 = "ident2", min.pct =min.pct, logfc.threshold = logfc.threshold)
GEMLI_tip[['DEG']]<-DEG

DEG_volcano_plot(GEMLI_tip, name1="sym_Tip", name2="asym_Tip")





# Extract endocrine cells

endocrine_cells <- GEMLI_items$cell_type$cell.ID[GEMLI_items$cell_type$cell.type %in% c("Early_trunk", "EP", "Alpha", "Beta")]


# Subset all relevant parts of the list using these cell IDs
GEMLI_endocrine <- list(
  gene_expression = GEMLI_items$gene_expression[, endocrine_cells, drop = FALSE],
  prediction = GEMLI_items$prediction[endocrine_cells, endocrine_cells, drop = FALSE],
  cell_type = GEMLI_items$cell_type[GEMLI_items$cell_type$cell.ID %in% endocrine_cells, ],
  cell_type_color = GEMLI_items$cell_type_color,
  predicted_lineages = GEMLI_items$predicted_lineages[endocrine_cells],
  predicted_lineage_table = GEMLI_items$predicted_lineage_table[GEMLI_items$predicted_lineage_table[, "cell.ID"] %in% endocrine_cells, ]
)

GEMLI_endocrine <- prediction_to_lineage_information(GEMLI_endocrine, cutoff=70)
cell_type_composition_plot(GEMLI_endocrine, ground_truth=F, cell_type_colors=T, type=c("upsetR")) 


# Extrace cell fate lineages
GEMLI_endocrine<-extract_cell_fate_lineages(GEMLI_endocrine, selection=c("EP", "Beta"), unique=FALSE, threshold = c(10,10))
GEMLI_endocrine[['cell_fate_analysis']][1:10,] 
table(GEMLI_endocrine[['cell_fate_analysis']]$cell.fate)

# Perform cell fate DEG calling
ident1 <- "sym_EP"
ident2 <- "asym_EP"
min.pct <- 0.05
logfc.threshold <- 0.1

GEMLI_Seurat<-CreateSeuratObject(GEMLI_endocrine[['gene_expression']], project = "SeuratProject", assay = "RNA")
Metadata<-GEMLI_endocrine[['cell_fate_analysis']];Metadata$ident<-NA
Metadata$ident[Metadata$cell.fate %in% ident1]<-"ident1"
Metadata$ident[Metadata$cell.fate%in%ident2]<-"ident2"
Meta<-as.data.frame(Metadata[,c(5)]); rownames(Meta)<-Metadata$cell.ID; colnames(Meta)<-c("cell.fate")
GEMLI_Seurat<-AddMetaData(GEMLI_Seurat, Meta, col.name = NULL); DefaultAssay(object = GEMLI_Seurat) <- "RNA"; Idents(GEMLI_Seurat) <- GEMLI_Seurat$cell.fate
GEMLI_Seurat <- NormalizeData(GEMLI_Seurat)
DEG <- FindMarkers(object = GEMLI_Seurat, ident.1 = "ident1", ident.2 = "ident2", min.pct =min.pct, logfc.threshold = logfc.threshold)
GEMLI_endocrine[['DEG']]<-DEG

DEG_volcano_plot(GEMLI_endocrine, name1="sym_EP", name2="asym_EP")
