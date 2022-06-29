library(qiime2R)
library(tidyverse)

# unpacking feature-table.qza file, extracting the ASV-table as data frame, and renaming appropriately
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
