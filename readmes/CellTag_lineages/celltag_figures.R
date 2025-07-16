# This script aims to plot the number of CellTag per cell and per sample

library(igraph)
library(proxy)
library(corrplot)
library(data.table)
library(CellTagR)
library(dplyr)
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

# Set the barcodes as row names
mat <- combined_data
rownames(mat) <- mat$Cell.BC
head(mat)

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



tags_per_day <- as.data.frame(tags_per_day)
tags_per_day$Day <- factor(rownames(tags_per_day), levels = c("iPSCs", "D2", "D6", "D8", "D10"))

ggplot(tags_per_day, aes(x = Day, y = tags_per_day)) +
  geom_col(fill = "steelblue") +
  labs(title = "Total Expressed CellTags per Day",
       x = "Day", y = "# CellTags Expressed") +
  theme_minimal()+
  theme(axis.text.x = element_text(hjust = 1, size = 10))




# count how many times each CellTag is expressed
# Step 1: Convert to long format
long_df <- mat %>%
  select(Day, 2:30) %>%
  pivot_longer(
    cols = -Day,
    names_to = "CellTag",
    values_to = "Count"
  )
long_df

# Step 2: Count how many cells have non-zero expression for each CellTag per day
summary_df <- long_df[long_df$Count > 0, ]
summary_df <- aggregate(Count ~ Day + CellTag, data = summary_df, FUN = length)
colnames(summary_df)[3] <- "CellCount"
summary_df

# Count non-zero CellTags per row (cell)
celltag_count_df <- data.frame(
  Cell.BC = mat[[1]],
  NumCellTags = rowSums(mat[, 2:(ncol(mat) - 2)] > 0)
)

# Remove cells with 0 expressed CellTags
celltag_count_df <- celltag_count_df[celltag_count_df$NumCellTags > 0, ]
celltag_count_df

ggplot(celltag_count_df, aes(x = reorder(Cell.BC, -NumCellTags), y = NumCellTags)) +
  geom_col(fill = "steelblue") +
  labs(title = "Number of Expressed CellTags per Cell",
       x = "Cell Barcode", y = "Number of CellTags Expressed") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))


ggplot(celltag_count_df, aes(y = NumCellTags)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  labs(title = "Distribution of Expressed CellTags per Cell",
       y = "# CellTags Expressed") +
  theme_minimal()




