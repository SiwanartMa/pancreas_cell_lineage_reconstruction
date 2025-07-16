# This script aims to make CellTag-Cell table for supplementary table S1

# Load required libraries
library(tidyverse)

# Load the csv file
df <- read.csv("/Users/mayongzhi/Desktop/researchProject/integration/metadata/celltag_counts.csv")
head(df)

# Keep only CellTags columns
df$Cell.BC <- NULL
df$n_tags <- NULL

# Rename the barcodes so that they align with the order during library prep
df <- df %>%
  mutate(X = sapply(strsplit(as.character(X), "\\+"), function(x) paste(c(x[2], x[1], x[-(1:2)]), collapse = "+")))

# Reshape long and filter for expressed CellTags
df_long <- df %>%
  pivot_longer(
    cols = -X,
    names_to = "CellTag",
    values_to = "count"
  ) %>%
  filter(count > 0)

# Reshape back — one row per CellTag, columns are barcodes
df_wide <- df_long %>%
  group_by(CellTag) %>%
  summarise(across(X, list)) %>%
  unnest_wider(X, names_sep = "_cell")

# View result
df_wide

# Export CSV
write.csv(df_wide, file = "~/Desktop/researchProject/integration/metadata/celltag_cells.csv")
