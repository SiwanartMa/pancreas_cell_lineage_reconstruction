#R
library(igraph)
library(proxy)
library(corrplot)
library(data.table)
library(CellTagR)
source("/Users/mayongzhi/Desktop/researchProject/integration/readmes/ScaleRNA_Rscript/CellTag source code/CellTagCloneCalling_Function.R")

# Vector of directories
dirs <- list.dirs(path = "~/Desktop/researchProject/CellTag_practice/CellTag_all/results/",
                  recursive = FALSE)

# Initialize list to store data
data_list <- list()

# Loop through each directory
for (dir in dirs) {
  # Extract sample ID from the directory name
  sample_id <- basename(dir)
  
  # Find the .Rds file (assumes only one matching RDS per directory)
  rds_file <- list.files(path = dir, pattern = "\\.celltag\\.matrix\\.Rds$", full.names = TRUE)
  
  # Skip if no RDS found
  if (length(rds_file) == 0) {
    message("No RDS file found in: ", dir)
    next
  }
  # Read and tag the data
  mat <- as.data.frame(readRDS(rds_file[1]))
  mat$Sample <- sample_id
  
  # Append to list
  data_list[[sample_id]] <- mat
}


# Get the union of all column names across all data frames
all_cols <- unique(unlist(lapply(data_list, names)))

# Ensure every data frame has all columns (fill missing with 0)
data_list_filled <- lapply(data_list, function(df) {
  missing_cols <- setdiff(all_cols, names(df))
  for (col in missing_cols) {
    df[[col]] <- 0
  }
  # Reorder columns to match
  df <- df[, all_cols]
  return(df)
})

# Combine all data into one data frame
combined_data <- do.call(rbind, data_list_filled)

# View result
dim(combined_data)
head(combined_data)

mat <- combined_data
mat <- mat[,-1]

# Extract differentiation day or iPSCs from the Sample column
mat$Day <- sub("CellTag-(D\\d+|iPSCs).*", "\\1", mat$Sample)
mat$Sample <- NULL

# Check the new column
table(mat$Day)
colnames(mat)

# Get only the CellTag columns (excluding "Day")
celltag_cols <- setdiff(colnames(mat), "Day")

# Create an empty vector to store results
tags_per_day <- setNames(numeric(0), character(0))

# Loop through each unique day
for (day in unique(mat$Day)) {
  # Subset cells for that day
  day_subset <- mat[mat$Day == day, celltag_cols]
  
  # Count how many CellTags are expressed (non-zero in any cell)
  expressed_tags <- colSums(day_subset > 0) > 0
  tags_per_day[day] <- sum(expressed_tags)
}

# View results
print(tags_per_day)

# Perform Single Cell Data Binarization
mat.bin$Day <- NULL
mat.bin <- SingleCellDataBinarization(celltag.dat = mat.bin, 1)
mat.bin

# Filter celltags by the whitelist
whitelist.path <- "/Users/mayongzhi/Desktop/researchProject/CellTag_practice/combined_celltags.csv"
mat.filt <- SingleCellDataWhitelist(celltag.dat = mat.bin,
                                    whitels.cell.tag.file = whitelist.path)

#mat.filt <- MetricBasedFiltering(whitelisted.celltag.data = mat.filt, cutoff = 0, comparison = "less")
mat.filt <- MetricBasedFiltering(whitelisted.celltag.data = mat.filt, cutoff = 1, comparison = "greater")

mat.sim <- JaccardAnalysis(whitelisted.celltag.data = mat.filt, plot.corr = TRUE, id = "all")

output.path <- "~/Desktop/researchProject/integration/outputs/ScaleRNA/celltag/"
mat.clones <- CloneCalling(Jaccard.Matrix = mat.sim, output.dir = output.path, output.filename = "all.clones.csv", correlation.cutoff = 0.1)

mat.clones

