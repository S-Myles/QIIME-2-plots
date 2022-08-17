# Loading packages
library(phyloseq)     # Data structure and functions for seq data
library(tidyverse)    # Data handling and all
library(microbiome)   # Contains CLR data transformation

### Importing data ###
# (See QIIME-2-plots/00-.../01-integrate_datasets_into_phyloseq.R for more info)
# Metadata
metadata <-  read_tsv("./data/metadata.tsv") %>%
  column_to_rownames(var = "sampleid")    # make sample ID the row names
# ASV table (Sample X Taxon, read counts)
ASV_table <- read_csv("./data/ASV_table.csv") %>%
  column_to_rownames('ASVs')      # Make ASV IDs row names
### Creating Phyloseq Object ###
ASV_table <- otu_table(ASV_table, taxa_are_rows = TRUE)         
metadata <- sample_data(metadata)                               
(physeq_data <-  merge_phyloseq(ASV_table, metadata)) 


### Plotting Shannon evenness boxplots ###
# There is lots to play with here to customize your viz. 
# Phyloseq plots are built atop of ggplot2, so most ggplot2 utilities can be used
plot_richness(physeq_data, x = "year", color="year", measure="Shannon") + 
  geom_boxplot(aes(group=year)) +
  theme(legend.position = "none") +
  coord_flip()
# Saving plot in /docs/
#ggsave("./docs/read_counts_Shannon.png")
