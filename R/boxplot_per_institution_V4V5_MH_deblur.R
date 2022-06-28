### Loading Packages and Setting Up Working Directory ###

# Run lines as needed to download packages you don't already have on your machine
#install.packages("BiocManager")
#BiocManager::install("phyloseq")
#or #install.packages("devtools")
#    devtools::install_github("joey711/phyloseq")
#install.packages("tidyverse")

# Loading packages
library(phyloseq)     # Data structure and functions for seq data
library(tidyverse)    # Data handling and all
library(microbiome)   # Contains CLR data transformation

setwd("C:/Users/18192/Documents/LaRoche_lab/GLOMICON/Inter_MetaBarcoding_Comparison")


### Pre-Processing All Data Sources ###

## Metadata
# import data and make sample ID the row names
metadata <-  read_tsv("./data/metadata.txt") %>%
  column_to_rownames(var = "sampleid")
# Specifying that DNA concentration should be a numeric value bc one 'N/A' renders it a character vector
metadata[['original DNA concentration (ng/uL)']] <- as.numeric(metadata[['original DNA concentration (ng/uL)']])
## ASV table (Sample X Taxon, read counts)
# Read ASV table but skip first row because all it is is '# Constructed from biom file'
ASV_table <- read_tsv("./data/V4V5_MH_deblur/feature-table.biom.tsv", skip=1) %>%
  as.data.frame() %>%             # Necessary format for phyloseq
  column_to_rownames('#OTU ID')   # Make ASV IDs row names


### Creating Phyloseq Object ###

# loading 18S V4 marker data, pre-processed with deblur, into respective phyloseq objects
ASV_table <- otu_table(ASV_table, taxa_are_rows = TRUE)         
metadata <- sample_data(metadata) 

# Merging all data into 1 global object
# (ASV table, Taxonomy table, and Metadata)
(physeq_data <-  merge_phyloseq(ASV_table, metadata))


### Plotting Chao1 boxplots for each providing institution
plot_richness(physeq_data, x = "provider", color="provider", measure="Chao1") + 
  geom_boxplot() +
  theme(legend.position = "none") +
  coord_flip()+ 
  ggtitle("V4V5 Amplicons processed with deblur (Chao1 index on read counts)")
# Sving plot in /doc/a_diversity
#ggsave("./doc/a_diversity/V4V5_deblur_Chao1.png")
### Plotting Shannon evenness boxplots for each providing institution
plot_richness(physeq_data, x = "provider", color="provider", measure="Shannon") + 
  geom_boxplot() +
  theme(legend.position = "none") +
  coord_flip() + 
  ggtitle("V4V5 Amplicons processed with deblur (Shannon index on read counts)")
# Sving plot in /doc/a_diversity
#ggsave("./doc/a_diversity/V4V5_deblur_Shannon.png")