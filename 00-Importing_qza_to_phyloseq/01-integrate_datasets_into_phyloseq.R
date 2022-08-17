# Loading packages
library(phyloseq)     # Data structure and functions for seq data
library(tidyverse)    # Data handling and all

# .qza files have already been transformed into .csv (See 00-.../00-...)

# Metadata
metadata <-  read_tsv("./data/metadata.tsv") %>%
  column_to_rownames(var = "sampleid")    # make sample ID the row names

# ASV table (Sample X Taxon, read counts)
ASV_table <- read_csv("./data/ASV_table.csv") %>%
  column_to_rownames('ASVs')      # Make ASV IDs row names

# Taxonomy reference table
taxonomy <- read_csv("./data/taxonomy.csv") %>%
  column_to_rownames('ASVs') %>%
  as.matrix()

### Creating Phyloseq Object ###

# loading data into respective phyloseq objects
ASV_table <- otu_table(ASV_table, taxa_are_rows = TRUE)         
metadata <- sample_data(metadata)  
taxonomy <- tax_table(taxonomy)

# Merging into 1 global object
# ASV IDs and sample names that don't match will get dropped
(physeq_data <-  merge_phyloseq(ASV_table, metadata, taxonomy)) 
# Check that your data was well loaded correctly in phyloseq.
# Encircling the above assignment between () will print the output assigned to physeq_data
