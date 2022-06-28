# QIIME-2-plots
A collection of convenience scripts to produce meaningful vizualisations from pre-processed metabarcoding data.
<br>

## Scope
The goal of this repository is to publish snippets of code commonly used in the downstream processing of amplicon data, often generated for microbiome and environmental DNA analysis. Following bioinformatics porcessing of the sequencing data into taxonomically assigned Amplicon Sequence Variants, few tools integrated in the processing software are flexible enough to accomodate the diverse needs of researchers. It's why many packages and libraries have been developped in higher-level programming languages to handle such tasks. My goal here is to bridge the gap between one common bioinformatics pipeline (QIIME) and many R packages to generate informative and flexible data vizualisations. 
<br>

### Amplicon Sequencing and pre-processing
This repository assumes that amplicon sequencing data has been pre-processed through a QIIME2 environment, as outlined in https://github.com/LangilleLab/microbiome_helper. Amplicons, processing tools employed within Qiime2, and parameter settings can vary. The output data expected (but not necessary for each applications here) are the following:
* Metadata.tsv
* ASV_table.qza
* Classification.qza

Example data is provided in the QIIME-2-Plots/data directory in the form of qiime2 artifacts (*.qza files). 


### QIIME-2-plots tools
