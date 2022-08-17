library(qiime2R)
library(tidyverse)

# unpacking classification.qza file, extracting the taxonomy data frame
taxonomy <- read_qza("data/classification.qza")['data'] %>% 
  data.frame() %>% 
  rename(ASVs = data.Feature.ID)  # Rename column with ASV IDs

# Cleanup taxonomic level prefixes (this will vary depending on your QIIME output)
# Here I will split the strings in data.Taxon column into 7 taxonomic levels (columns)

# What I'm starting with: 'd__Eukaryota; p__Diatomea; c__Mediophyceae; o__Mediophyceae; f__Mediophyceae; g__Thalassiosira'
# What I want:| Eukaryota | Diatomea | Mediophyceae | Mediophyceae | Mediophyceae | Thalassiosira |         | 

# remove prefixes
taxonomy$data.Taxon <- gsub('.__', '', taxonomy$data.Taxon)
# Make a list of the  maximum taxonomic levels available in data.Taxon column (an artifact of your classifier)
tax_levels <- c('domain', 'phyla', 'class', 'order', 'family', 'genus', 'specie')
# Split 'data.taxon' column into 7 taxonomic levels (columns)
taxonomy[tax_levels] <- str_split_fixed(taxonomy$data.Taxon, '; ', 7)

# Tidying dataset (I only care for the ASV IDs and taxonomic level columns)
taxonomy <- taxonomy[c('ASVs', tax_levels)]

# writing .csv file into 'data/'
write.csv(taxonomy, 'data/taxonomy.csv', row.names = F)