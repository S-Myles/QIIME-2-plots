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

# This data set is split into two pools of eukaryotic and prokaryotic sequences. So the following methods will be applied in parallel where needed.


### Pre-Processing All Data Sources ###

## Metadata
# import data and make sample ID the row names
metadata <-  read_tsv("./data/metadata.txt") %>%
  column_to_rownames(var = "sampleid")
# Specifying that DNA concentration should be a numeric value bc one 'N/A' renders it a character vector
metadata[['original DNA concentration (ng/uL)']] <- as.numeric(metadata[['original DNA concentration (ng/uL)']])
## ASV table (Sample X Taxon, read counts)
# Read ASV table but skip first row because all it is is '# Constructed from biom file'
# PROKs
PROK_ASV_table <- read_tsv("./data/V4V5_BB_DADA2/PROK/feature-table.biom.tsv", skip=1) %>%
  as.data.frame() %>%             # Necessary format for phyloseq
  column_to_rownames('#OTU ID')   # Make ASV IDs row names
# EUKs
EUK_ASV_table <- read_tsv("./data/V4V5_BB_DADA2/EUK/feature-table.biom.tsv", skip=1) %>%
  as.data.frame() %>%             # Necessary format for phyloseq
  column_to_rownames('#OTU ID')   # Make ASV IDs row names

### Creating Phyloseq Object ###
## loading universal V4V5 marker data, pre-processed with DADA2, into respective phyloseq objects
# PROKs
PROK_ASV_table <- otu_table(PROK_ASV_table, taxa_are_rows = TRUE)         
# EUKs
EUK_ASV_table <- otu_table(EUK_ASV_table, taxa_are_rows = TRUE)         
# Metadata is the same for all
metadata <- sample_data(metadata)  
## Merging all data into 1 global object
#(ASV table, Taxonomy table, and Metadata)
# PROKs
(PROK_V4V5_BB_DADA2_data <-  merge_phyloseq(PROK_ASV_table, metadata))
# EUKs
(EUK_V4V5_BB_DADA2_data <-  merge_phyloseq(EUK_ASV_table, metadata))


### Plotting Chao1 boxplots for each providing institution
# PROKs
plot_richness(PROK_V4V5_BB_DADA2_data, x = "provider", color="provider", measure="Chao1") + 
  geom_boxplot() +
  theme(legend.position = "none") +
  coord_flip()+ 
  ggtitle("V4V5 16S Amplicons processed with Furhman's DADA2 (Chao1 index on read counts)")
# Sving plot in /doc/a_diversity
#ggsave("./doc/a_diversity/V4V5_BB_DADA2_PROKs_Chao1.png")

# EUKs
plot_richness(EUK_V4V5_BB_DADA2_data, x = "provider", color="provider", measure="Chao1") + 
  geom_boxplot() +
  theme(legend.position = "none") +
  coord_flip()+ 
  ggtitle("V4V5 18S Amplicons processed with Furhman's DADA2 (Chao1 index on read counts)")
# Sving plot in /doc/a_diversity
#ggsave("./doc/a_diversity/V4V5_BB_DADA2_EUKs_Chao1.png")

### Plotting Shannon evenness boxplots for each providing institution
# PROKs
plot_richness(PROK_V4V5_BB_DADA2_data, x = "provider", color="provider", measure="Shannon") + 
  geom_boxplot() +
  theme(legend.position = "none") +
  coord_flip() + 
  ggtitle("V4V5 16S Amplicons processed with Furhman's DADA2 (Shannon index on read counts)")
# Sving plot in /doc/a_diversity
#ggsave("./doc/a_diversity/V4V5_BB_DADA2_PROKs_Shannon.png")

# EUKs
plot_richness(EUK_V4V5_BB_DADA2_data, x = "provider", color="provider", measure="Shannon") + 
  geom_boxplot() +
  theme(legend.position = "none") +
  coord_flip() + 
  ggtitle("V4V5 18S Amplicons processed with Furhman's DADA2 (Shannon index on read counts)")
# Sving plot in /doc/a_diversity
#ggsave("./doc/a_diversity/V4V5_BB_DADA2_EUKs_Shannon.png")
