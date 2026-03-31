# RNAseq-DCM-GeneExpression-Analysis
RNA-seq analysis of dilated cardiomyopathy (DCM) vs control heart samples using GEO dataset GSE141910. Includes PCA, UMAP, volcano plot, and heatmap visualizations in R.

# Introduction
DCM is a condition where the heart becomes enlarged and weakened, making it harder to pump blood effectively. As the left ventricle stretches, the heart’s ability to supply the body with oxygen decreases, often leading to heart failure. It is one of the leading causes of systolic heart failure and can affect people at relatively young ages. Once symptoms appear, the disease can be serious, with significant risk of death within a few years. DCM also increases the risk of dangerous heart rhythms that can cause sudden cardiac death. Because of its severity and the fact that its causes are often complex or unclear, it is an important condition to study in order to improve treatment and outcomes.

# Research Question
Which genes and molecular pathways differ between males and females with dilated cardiomyopathy (DCM) compared to non-failing hearts?

# Hypothesis
DCM hearts exhibit distinct gene expression and pathway profiles compared to non-failing hearts, characterized by increased inflammatory signaling and activation of cardiac remodeling programs, with these effects being more pronounced in males than in females. 

# RNA-seq Analysis of Dilated Cardiomyopathy (DCM)

This project analyzes RNA-seq gene expression data from the GEO dataset **GSE141910** to identify differential gene expression between **Dilated Cardiomyopathy (DCM)** and **Non-failing human hearts**. Furthermore, this project will evluate if these transcriptomic differences are **sex-specific**.  

# Dataset: RNA sequencing of the left ventricle from non-failing donors and heart failure samples from the MAGNet consortium

Human left ventricular RNA-seq data. Includes non-failing hearts and multiple cardiomyopathy subtypes (including DCM)
n = 305
- 152 patients with non-failing donors
    - 70 Males
    - 82 Females
- 153 patients with dilated cardiomyopathy
    - 92 Males
    - 61 Females
- Overall sex distribution:
    - 53% Male / 47% Female


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
