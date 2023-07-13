library(tidyverse)

# Import data
ASV_table <- read.csv("feature-table.tsv", sep="\t", row.names=1)

# If ASV table is oriented with samples as columns and ASVs as rows,

# ASV Abundance filtering
ASV_table_filt <- ASV_table[rowSums(ASV_table)>=100,] 

# ASV Prevalence filtering
rows_to_keep <- rowSums(ASV_table_filt != 0) >= 5 #Find rows (ASVs) that have more than 5 occurrences
ASV_table_filt <- ASV_table_filt[rows_to_keep,] # Keep only those rows



######################################################################
# If ASV table is transposed with ASVs as columns and samples as rows,
######################################################################
# If you want to transpose yourself,
#ASV_table_transposed <- t(ASV_table) %>% as_tibble(rownames= NA)

# ASV Abundance filtering
#ASV_table_filt <- ASV_table[colSums(ASV_table)>=100]

# ASV Prevalence filtering
#columns_to_keep <- colSums(ASV_table_filt != 0) >= 5 #Find columns (ASVs) that have more than 5 occurrences
#ASV_table_filt <- ASV_table_filt[, columns_to_keep] # Keep only those columns
