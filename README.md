# QIIME-2-plots
A collection of R scripts to handle pre-processed metabarcoding data and produce meaningful vizualisations.
<br>

## Scope
The goal of this repository is to publish snippets of open source code that I regularly use for downstream processing of amplicon sequencing data. This type of data is often generated for microbiome and environmental DNA analysis. Following the processing of raw sequencing data through a typical bioinformatics pipeline ([QIIME2](https://qiime2.org/)), few of the integrated tools are flexible enough to accomodate the diverse research needs. It's why many packages and libraries have been developped and published in higher-level programming languages to handle such tasks. 
<br>
My goal here is to bridge the gap between one bioinformatics pipeline and many R packages to generate informative and flexible data vizualisations. 
<br>

## Data
### Amplicon Sequencing and pre-processing
This repository assumes that amplicon sequencing data generated from an Illumina MiSeq platform has been pre-processed with the QIIME2 environment, as outlined in https://github.com/LangilleLab/microbiome_helper. Amplicon sources, processing tools employed within Qiime2, and parameter settings along the pre-processing pipeline will have an impact on results and can vary. 
<br>

### Input data for QIIME-2-plots
The data expected (not all sources are used in all functions and additional inputs could easily be added):
* Metadata.tsv
> File imported into QIIME2 and used during pre-processing. <br>
> First column is 'sampleid' (or a similar title accepted by the QIIME2 format) and must identically match columns of the feature (ASV) table. <br>
* ASV_table.qza
> 1st column 'Feature-ID' are the ASV IDs and must match ASV IDs of classification.qza. <br>
> Next columns are samples with read counts for each ASV detected. <br>
> Samples whose string must match strings of column 'sampleid' in Metadata.tsv. <br>
* Classification.qza
> 1st column 'Feature-ID' are the ASV IDs and must match ASV IDs of ASV_table.qza. <br> 
> 2nd column will be the taxonomic assignement which will vary depending on the classifier used within the QIIME2 environment. <br>

Input data used to create these scripts have been removed from this repository for confidentiality. When inputing your data, you will need to adapt scripts to your character patterns and purposes. 

`*` Checking your data on the console as lines of code get run will diagnose issues and save you a headache. 
<br>

## Tools
