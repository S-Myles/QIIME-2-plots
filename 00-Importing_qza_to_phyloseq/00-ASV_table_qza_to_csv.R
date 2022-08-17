library(qiime2R)
library(tidyverse)

# unpacking feature-table.qza file and extracting the ASV-table as data frame
# You will likely need to rename the columns of this table (samples) to match metadata file.
# Your patterns will deviate from this, and in need, lookup 'RegEx patterns' for guidance.
ASV_table <- read_qza("data/dada2_table_final.qza")['data'] %>% 
  data.frame() %>% 
  rename_with(                   # Renaming column names to match metadata file
    function(x){
    x <- gsub('\\.', '-', x)     # Changing all '.' by '-'
    gsub('data-', '', x)         # removing all instances of 'data-'
      })

ASV_table <- cbind(ASVs=rownames(ASV_table), ASV_table) 

# writing .csv file into 'data/'
write.csv(ASV_table, 'data/ASV_table.csv', row.names = F)
