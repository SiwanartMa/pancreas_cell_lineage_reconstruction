# Install packages
library(Seurat)
library(tidyverse)
library(patchwork)

# Load the 3 Seurat objects
# Increase the memory
mem.maxVSize(vsize = 90000)

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

# Plot the time point each SCT assay covers
library(ggplot2)
library(dplyr)
library(tidyr)

# Create a long-format dataframe
cluster_cols <- c("SCT_snn_res.0.2", "SCT_snn_res.0.1", "SCT_snn_res.0.05")
meta_data <- human12_20W@meta.data

long_df <- meta_data %>%
  select(Gestation_Age, all_of(cluster_cols)) %>%
  pivot_longer(cols = all_of(cluster_cols), names_to = "resolution", values_to = "cluster") %>%
  mutate(cluster = as.factor(cluster), resolution = as.factor(resolution), gestation_age = as.factor(Gestation_Age))

# Plot
ggplot(long_df, aes(x = cluster, fill = gestation_age)) +
  geom_bar(position = "fill") +
  facet_wrap(~ resolution, scales = "free_x") +
  labs(x = "Cluster", y = "Proportion", fill = "Gestation Age") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("../outputs/multiSCT_cover.png")

# Remove SCT_snn_res.0.1 and SCT_snn_res.0.05
human12_20W@meta.data <- human12_20W@meta.data %>%
  select(-SCT_snn_res.0.1, -SCT_snn_res.0.05)

# Set the SCT assay as the default assay
DefaultAssay(human12_20W) <- "SCT"

# Filter the object based on the criteria in the paper
human12_20W<- subset(human12_20W, subset = nCount_RNA > 500 & nFeature_RNA > 250 & log10GenesPerUMI > 0.75 & percent.mt < 25)

# Combine data in 13PCW and 13PCW_1
human12_20W@meta.data$Gestation_Age[human12_20W@meta.data$Gestation_Age == "13PCW_1"] <- "13PCW"

# Add time column in each object
human4_6W$time <- "4_6"
human7_11W$time <- "7_11"
human12_20W$time <- "12_20"


# Modify the "Gestation_Age" column in metadata
human12_20W@meta.data$Gestation_Age <- human12_20W@meta.data$Gestation_Age %>%
  as.character() %>%  # Ensure it's a character vector
  gsub("PCW", "", .) %>%  # Remove "PCW"
  paste0("W", .)  # Add "W" at the beginning

# Rename "Gestation_Age" with "days"
colnames(human12_20W@meta.data)[colnames(human12_20W@meta.data) == "Gestation_Age"] <- "days"

# Make a list storing all Seurat objects
obj_list <- list(human4_6W, human7_11W, human12_20W)

# Merge
combined <- merge(human4_6W, y = list(human7_11W, human12_20W))
DefaultAssay(combined) <- "SCT"
saveRDS(combined, file = "/Users/mayongzhi/Desktop/researchProject/integration/outputs/merged_SCT.rds")

# Make sure each object has only one SCT model
human4_6W@assays$SCT #1assay
human7_11W@assays$SCT #1assay
human12_20W@assays$SCT #4assays

# human12_20W has 4 SCTassays
# Set it to NULL and rerun it
DefaultAssay(human12_20W) <- "RNA"
human12_20W[["SCT"]] <- NULL 
options(future.globals.maxSize = 8 * 1024^3) # Set memory limit to 8 GB
human12_20W <- SCTransform(human12_20W, verbose = TRUE) 
human12_20W@assays$SCT@SCTModel.list # 1assay
DefaultAssay(human12_20W) <- "SCT"

# Set variable features and run PCA
human12_20W <- FindVariableFeatures(human12_20W)

# Make the object list again
obj_list <- list(human4_6W, human7_11W, human12_20W)

# Select integration features
features <- SelectIntegrationFeatures(
  object.list = obj_list, nfeatures = 3000,
  verbose = TRUE)

VariableFeatures(combined) <- features
combined <- RunPCA(combined, assay = "SCT", features = features)

# Make elbow plot
ElbowPlot(combined)
#ggsave("/Users/mayongzhi/Desktop/researchProject/integration/outputs/unintegrated_figure/elbow_combined.png")

combined <- IntegrateLayers(
  object = combined, 
  method = CCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca",
  new.reduction = "integrated.pca",
  verbose = TRUE
)
