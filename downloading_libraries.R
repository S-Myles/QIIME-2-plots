# Run lines as needed to download packages you don't already have on your machine

install.packages("tidyverse")

# With BiocManager
install.packages("BiocManager")

BiocManager::install("phyloseq")
#or #install.packages("devtools")
#    devtools::install_github("joey711/phyloseq")
BiocManager::install("microbiome")

# With remotes
install.packages("remotes")

remotes::install_github("jbisanz/qiime2R")