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

############ have a look at the read distribution
#Make a data frame with a column for the read counts of each sample##### 
sample_sum_df <- data.frame(sum = sample_sums(physeq_data))
# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 5000) +
  ggtitle("Distribution of sequencing depths of all samples") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())


# Sample read depth boxplots
ggplot(sample_sum_df,aes(x="", y=sum)) + geom_violin(trim=FALSE) + 
  geom_boxplot(width=0.25) +
  labs(
    x = "Samples",
    y = "ASV Read Depth"
  )