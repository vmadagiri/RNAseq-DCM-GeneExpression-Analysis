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

![PCA](PCA Plot - Sex Comparison.png)
<img width="1552" height="1222" alt="PCA Plot - Sex Comparison" src="https://github.com/user-attachments/assets/57057b5e-1f44-4a9c-82cc-851477eac122" />

### UMAP Plot

![UMAP](UMAP Plot - Sex Comparison.png)
<img width="1552" height="1222" alt="UMAP Plot - Sex Comparison" src="https://github.com/user-attachments/assets/6d7b422a-15db-481b-9550-0d09f0ed4262" />

### Volcano Plot (Female vs Male DCM)

![Volcano](Volcano Plot - Female vs Male DCM.png)
<img width="1552" height="1222" alt="Volcano Plot - Female vs Male DCM" src="https://github.com/user-attachments/assets/187b8d0e-1e44-4dbb-bd6b-251c8d770b3c" />

### Volcano Plot (Female vs Male Control)

![Volcano](Volcano Plot - Female vs Male control.png)<img width="1552" height="1222" alt="Rplot" src="https://github.com/user-attachments/assets/2b3385cd-bf17-4751-ad6c-1435dd842458" />


![Uploading Rplot.png…]()

### PC1 Distribution

![PC1](PC1 Distribution by Group.png)
<img width="1552" height="1222" alt="PC1 Distribution by Group" src="https://github.com/user-attachments/assets/e27a32a3-c3c7-42a1-961b-8663568651e3" />

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

## Discussion

The results from this analysis show that gene expression patterns in DCM differ from non-failing hearts, although the separation is not completely clear. In the PCA plot, there is some overlap between groups, which suggests that disease status is not the only factor influencing gene expression. This makes sense given how complex DCM is and how much variability exists between patients.

The UMAP plot shows slightly clearer clustering compared to PCA, which suggests that there may be nonlinear patterns in the data that PCA is not capturing as well. Even so, the groups are still not perfectly separated, reinforcing the idea that multiple factors are contributing to the observed variation.

The differential expression results highlight genes that are upregulated and downregulated in DCM. Many of these genes are likely involved in processes such as inflammation and cardiac remodeling, which are commonly associated with heart failure. This aligns with what is already known about the disease.

There also appear to be differences between male and female samples, although these differences are not extremely distinct. This suggests that sex may play a role in how DCM develops or progresses, but it is probably not the only factor. More focused analysis would be needed to better understand these differences.

Overall, the results suggest that DCM is associated with measurable changes in gene expression, but the patterns are complex and not driven by a single variable. This highlights the importance of considering multiple biological and clinical factors when studying the disease.

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

## Limitations

* This analysis uses bulk RNA-seq data, so it doesn’t tell us which specific cell types the gene expression changes are coming from
*  There is a lot of variability between patients (like disease stage or treatment history), and that isn’t fully accounted for here
*  The sample sizes between groups aren’t perfectly balanced, which could affect some of the comparisons
*  All results are based on a single dataset, so they would need to be validated using additional datasets to confirm the findings  

---

## Conclusion

This study demonstrates that DCM is associated with meaningful transcriptomic changes, particularly in pathways related to inflammation and remodeling. The observed sex-specific differences suggest that personalized approaches may be important in understanding and treating the disease.

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
