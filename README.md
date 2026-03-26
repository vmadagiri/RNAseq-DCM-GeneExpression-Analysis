# RNAseq-DCM-GeneExpression-Analysis
RNA-seq analysis of dilated cardiomyopathy vs control heart samples using GEO dataset GSE141910. Includes PCA, UMAP, volcano plot, and heatmap visualizations in R.

# Introduction
Explain what DCM is and why it is important to study. Talk about mortality rates etc. 

# Research Question
Which genes and molecular pathways differ between males and females with dilated cardiomyopathy (DCM) compared to non-failing hearts?

# Hypothesis

# RNA-seq Analysis of Dilated Cardiomyopathy (DCM)

This project analyzes RNA-seq gene expression data from the GEO dataset **GSE141910** to identify differential gene expression between **Dilated Cardiomyopathy (DCM)** and **Non-failing human hearts**. Furthermore, this project will evluate if these transcriptomic differences are **sex-specific**.  

# Dataset


## Methods

The analysis pipeline includes:

- Data filtering of low-count genes
- Log transformation of counts
- Principal Component Analysis (PCA)
- UMAP dimensionality reduction
- Differential gene expression analysis using **limma**
- Volcano plot visualization
- Heatmap of top differentially expressed genes

## Dataset

GEO Accession:  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141910

## Visualizations

### PCA Plot
![PCA](PCA_GSE141910.png)

### UMAP Plot
![UMAP](UMAP_GSE141910.png)

### Volcano Plot
![Volcano](Volcano_GSE141910.png)

### Heatmap
![Heatmap](Heatmap_GSE141910.png)

## Tools Used

- R
- GEOquery
- limma
- umap
- data.table

## Author

Vaishnavi Madagiri  
Bioinformatics Major — Virginia Commonwealth University

Ammar Mohiuddin
Biology and Bioinformatics Major - Virginia Commonwealth University 
