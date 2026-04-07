# RNAseq-DCM-GeneExpression-Analysis

### RNA-seq Analysis using GSE141910

Exploring gene expression differences in dilated cardiomyopathy (DCM) with a focus on sex-specific effects

---

## Introduction

Dilated cardiomyopathy (DCM) is a condition where the heart becomes enlarged and weakened, reducing its ability to pump blood effectively. This can lead to heart failure, arrhythmias, and in severe cases, sudden cardiac death. DCM is one of the leading causes of heart transplantation and can affect individuals at relatively young ages.

The causes of DCM are complex and include genetic mutations, environmental factors, and unknown contributors. Because of this, studying the disease at the molecular level is important for understanding what is actually happening inside the heart.

This project uses RNA-seq data to analyze gene expression changes in DCM compared to non-failing hearts. In addition, this project specifically explores whether these changes differ between males and females.

---

## Research Question

Which genes and molecular pathways differ between males and females with dilated cardiomyopathy (DCM) compared to non-failing hearts?

---

## Hypothesis

DCM hearts show distinct gene expression and pathway changes compared to non-failing hearts, with stronger inflammatory and remodeling signals. These effects are expected to differ between males and females.

---

## Dataset: GSE141910

RNA sequencing of the left ventricle from non-failing donors and heart failure samples from the MAGNet consortium.

* Total samples: 305
* Non-failing: 152

  * 70 males
  * 82 females
* DCM: 153

  * 92 males
  * 61 females

GEO Accession:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141910

---

## Methods

The analysis pipeline includes:

* Filtering low-count genes
* Log transformation of counts
* Principal Component Analysis (PCA)
* UMAP clustering
* Differential gene expression using limma
* Volcano plot visualization
* Sex-based comparison of DCM samples
* KEGG pathway enrichment analysis

---

## Visualizations

### PCA Plot

![PCA](PCAplot.png)

### UMAP Plot

![UMAP](UMAP Plot - Sex Comparison.png)

### Volcano Plot (Female vs Male DCM)

![Volcano](Volcano Plot - Female vs Male DCM.png)

### PC1 Distribution

![PC1](PC1 Distribution by Group.png)

---

## 🧬 KEGG Pathway Analysis

KEGG enrichment analysis was performed on significantly differentially expressed genes to understand biological pathways involved in sex differences in DCM.

### KEGG Summary Barplot

![KEGG](sex_differences_analysis/KEGG_Summary_Barplot.png)

*Barplot showing enriched pathways related to inflammation and metabolic processes in male vs female DCM samples.*

---

## 🔬 Key Findings

* DCM samples separate clearly from controls in PCA and UMAP
* There are noticeable gene expression differences between male and female DCM hearts
* Several genes are significantly upregulated in DCM
* KEGG analysis highlights pathways related to inflammation and metabolism
* These results suggest sex-specific biological differences in DCM progression

---

## Tools Used

* R
* GEOquery
* limma
* clusterProfiler
* org.Hs.eg.db
* umap
* data.table

---

## Author

Vaishnavi Madagiri
Bioinformatics Major — Virginia Commonwealth University

Ammar Mohiuddin
Biology and Bioinformatics — Virginia Commonwealth University

Harrish Ganesh
Biology and Bioinformatics — Virginia Commonwealth University

Haneia Nemati
Bioinformatics — Virginia Commonwealth University
