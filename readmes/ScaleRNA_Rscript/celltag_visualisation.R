library(tidyverse)
library(foreach)
library(networkD3)

source("/Users/mayongzhi/Desktop/researchProject/integration/readmes/ScaleRNA_Rscript/CellTag source code/function_source_for_network_construction.R")
source("/Users/mayongzhi/Desktop/researchProject/integration/readmes/ScaleRNA_Rscript/CellTag source code/function_source_for_network_visualization.R")

clones <- read_csv("~/Desktop/researchProject/integration/outputs/ScaleRNA/celltag/all.clones.csv")
colnames(clones)[1] <- "CellTagV1"
clones
clone.cells <- unique(c(clones$cell.barcode))
celltag_data <- data.frame(clone.cells, row.names = clone.cells)
celltag_data
celltag_data$CellTagV1 <- NA
celltag_data[clones$cell.barcode, "CellTagV1"] <- clones$CellTagV1
#celltag_data <- celltag_data[, -1]
row.names(celltag_data) <- paste0(rownames(celltag_data), "-1")

celltag_data $ clone.cells <- NULL
celltag_data$CellTagV2 <- NA
celltag_data$CellTagV3 <- NA
colnames(celltag_data) <- c("CellTagV1", "CellTagV2", "CellTagV3")
celltag_data

linkList <- convertCellTagMatrix2LinkList(celltag_data)
Nodes <- getNodesfromLinkList(linkList)

drawSubnet(celltag_data, tag = "CellTagV1", linkList = linkList, Nodes = Nodes)

bar.data <- celltag_data
bar.data$Cell.BC <- row.names(bar.data)
bar.data <- gather(bar.data, key = "CellTag", value = "Clone", 1:3, na.rm = FALSE)
ggplot(data = bar.data) + geom_bar(mapping = aes(x = CellTag, fill = factor(Clone)), position = "fill", show.legend = FALSE) + scale_y_continuous(labels = scales::percent_format())
