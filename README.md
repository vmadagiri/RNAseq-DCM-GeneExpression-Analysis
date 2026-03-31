# RNAseq-DCM-GeneExpression-Analysis
RNA-seq analysis of dilated cardiomyopathy (DCM) vs control heart samples using GEO dataset GSE141910. Includes PCA, UMAP, volcano plot, and heatmap visualizations in R.

# Introduction
DCM is a condition where the heart becomes enlarged and weakened, making it harder to pump blood effectively. As the left ventricle stretches, the heart’s ability to circulate oxygen-rich blood decreases, often leading to heart failure. DCM is one of the leading causes of systolic heart failure and can affect individuals at relatively young ages, not just older adults. Once symptoms develop, the disease can be serious, with a significant risk of mortality within a few years depending on severity and access to treatment. Patients are also at risk for dangerous arrhythmias, which can lead to sudden cardiac death. Because of its severity and unpredictable progression, DCM remains an important focus of medical research.

The causes of DCM are complex and often multifactorial, including genetic mutations, environmental exposures, and underlying health conditions. In many cases, the exact cause is still unknown, which makes it difficult to predict disease progression or tailor treatments. This highlights the need for research that looks beyond clinical symptoms and focuses on the underlying biological mechanisms of the disease.

In this project, we focus on understanding DCM at the molecular level using RNA sequencing (RNA-seq) data. Rather than only observing clinical outcomes, this approach allows us to analyze gene expression patterns in heart tissue and identify which genes and pathways are altered in DCM compared to non-failing hearts. By comparing these differences, we can gain insight into the biological processes driving the disease. Overall, this project aims to contribute to a better understanding of DCM and support future efforts toward more targeted and personalized treatments.
# Research Question
Which genes and molecular pathways differ between males and females with dilated cardiomyopathy (DCM) compared to non-failing hearts?

# Hypothesis
DCM hearts exhibit distinct gene expression and pathway profiles compared to non-failing hearts, characterized by increased inflammatory signaling and activation of cardiac remodeling programs, with these effects being more pronounced in males than in females. 

# RNA-seq Analysis of Dilated Cardiomyopathy (DCM)

This project analyzes RNA-seq gene expression data from the GEO dataset **GSE141910** to identify differential gene expression between **Dilated Cardiomyopathy (DCM)** and **Non-failing human hearts**. Furthermore, this project will evluate if these transcriptomic differences are **sex-specific**.  

# Dataset: GSE141910

**RNA sequencing of the left ventricle from non-failing donors and heart failure samples from the MAGNet consortium**

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
 
**GEO Accession:**
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141910


## Methods

The analysis pipeline includes:

- Data filtering of low-count genes
- Log transformation of counts
- Principal Component Analysis (PCA)
- UMAP dimensionality reduction
- Differential gene expression analysis using **limma**
- Volcano plot visualization
- Heatmap of top differentially expressed genes

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
