library(tidyverse)
library(phyloseq)

# Metadata
metadata <-  read_tsv("data/metadata.tsv") %>%
  column_to_rownames(var = "sampleid")    # make sample ID the row names

# ASV table (Sample X Taxon, read counts)
ASV_table <- read_csv("data/ASV_table.csv") %>%
  column_to_rownames('ASVs')      # Make ASV IDs row names

### Creating Phyloseq Object ###
# loading data into respective phyloseq objects
ASV_table <- otu_table(ASV_table, taxa_are_rows = TRUE)         
metadata <- sample_data(metadata)                               
# Merging into 1 global object
(physeq_data <-  merge_phyloseq(ASV_table, metadata)) 

# compute ASV read counts of each sample into a data frame 
sample_sum_df <- data.frame(sum = sample_sums(physeq_data))

# Find outliers
out <- boxplot.stats(sample_sum_df$sum)$out
out_ind <- which(sample_sum_df$sum %in% c(out))
outliers <- row.names(metadata[out_ind])

# remove samples with outlier read depths
(physeq_data_no_out <- prune_samples(!(sample_names(physeq_data) %in% outliers), physeq_data))

# remove samples with arbitrary read depths threshold
(physeq_data_no_out <- prune_samples(sample_sums(physeq_data_no_out) >= 7000 , physeq_data_no_out))
