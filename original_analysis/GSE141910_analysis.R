# install.packages("BiocManager")
# install.packages("car")
# install.packages("umap")
# BiocManager::install("GEOquery")
# BiocManager::install("Biobase")
# BiocManager::install("limma")

library(data.table)
library(GEOquery)
library(Biobase)
library(limma)
library(umap)
library(car)

# -------------------------
# Load dataset
# -------------------------

urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(
  urld,
  "acc=GSE141910",
  "file=GSE141910_raw_counts_GRCh38.p13_NCBI.tsv.gz",
  sep = "&"
)

tbl <- as.matrix(
  data.table::fread(path, header = TRUE, colClasses = "integer"),
  rownames = 1
)

dim(tbl)
head(tbl[, 1:5])

# -------------------------
# Filter low-count genes
# -------------------------

keep <- rowSums(tbl >= 10) >= 2
tbl <- tbl[keep, ]

# -------------------------
# Log transform
# -------------------------

dat <- log10(tbl + 1)

# -------------------------
# Group labels
# -------------------------

group <- factor(c(rep("DCM", 166), rep("Control", 190)))
group_colors <- ifelse(group == "DCM", "red", "blue")

# -------------------------
# Boxplot
# -------------------------

boxplot(
  dat,
  boxwex = 0.5,
  notch = FALSE,
  outline = FALSE,
  las = 2,
  cex.axis = 0.5,
  main = "GSE141910 RNA-seq Distribution",
  ylab = "log10(count + 1)"
)

# -------------------------
# UMAP
# -------------------------

dat2 <- dat[!duplicated(dat), ]
ump <- umap::umap(t(dat2), n_neighbors = 15, random_state = 123)

plot(
  ump$layout,
  main = "GSE141910 UMAP Plot",
  xlab = "UMAP1",
  ylab = "UMAP2",
  pch = 20,
  cex = 1.2,
  col = group_colors
)

legend(
  "topright",
  legend = c("DCM", "Control"),
  col = c("red", "blue"),
  pch = 20
)

# -------------------------
# PCA Plot
# -------------------------

pca <- prcomp(t(dat), scale. = TRUE)

plot(
  pca$x[,1],
  pca$x[,2],
  col = group_colors,
  pch = 19,
  xlab = "PC1",
  ylab = "PC2",
  main = "PCA Plot - DCM vs Control"
)

legend(
  "topright",
  legend = c("DCM", "Control"),
  col = c("red", "blue"),
  pch = 19
)

# -------------------------
# Differential Expression
# -------------------------

design <- model.matrix(~group)

fit <- lmFit(dat, design)
fit <- eBayes(fit)

topGenes <- topTable(fit, coef = 2, number = Inf)

head(topGenes)

# -------------------------
# Volcano Plot
# -------------------------

plot(
  topGenes$logFC,
  -log10(topGenes$P.Value),
  pch = 20,
  col = ifelse(topGenes$adj.P.Val < 0.05, "red", "grey"),
  main = "Volcano Plot",
  xlab = "log2 Fold Change",
  ylab = "-log10 P-value"
)

abline(h = -log10(0.05), col = "blue", lty = 2)
abline(v = c(-0.1, 0.1), col = "darkgreen", lty = 2)

top10 <- head(topGenes[order(topGenes$adj.P.Val), ], 10)

text(
  top10$logFC,
  -log10(top10$P.Value),
  labels = rownames(top10),
  pos = 3,
  cex = 0.6
)

# -------------------------
# Heatmap of top genes
# -------------------------

top <- rownames(topGenes)[1:30]
sample_colors <- ifelse(group == "DCM", "red", "blue")

heatmap(
  dat[top, ],
  Colv = NA,
  scale = "row",
  labCol = FALSE,
  ColSideColors = sample_colors,
  col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  main = "Top Differentially Expressed Genes"
)

# -------------------------
# Save plots as PNG files
# -------------------------

png("Boxplot_GSE141910.png", width = 1000, height = 700)
boxplot(
  dat,
  boxwex = 0.5,
  notch = FALSE,
  outline = FALSE,
  las = 2,
  cex.axis = 0.5,
  main = "GSE141910 RNA-seq Distribution",
  ylab = "log10(count + 1)"
)
dev.off()

png("UMAP_GSE141910.png", width = 900, height = 700)
plot(
  ump$layout,
  main = "GSE141910 UMAP Plot",
  xlab = "UMAP1",
  ylab = "UMAP2",
  pch = 20,
  cex = 1.2,
  col = group_colors
)
legend(
  "topright",
  legend = c("DCM", "Control"),
  col = c("red", "blue"),
  pch = 20
)
dev.off()

png("PCA_GSE141910.png", width = 900, height = 700)
plot(
  pca$x[,1],
  pca$x[,2],
  col = group_colors,
  pch = 19,
  xlab = "PC1",
  ylab = "PC2",
  main = "PCA Plot - DCM vs Control"
)
legend(
  "topright",
  legend = c("DCM", "Control"),
  col = c("red", "blue"),
  pch = 19
)
dev.off()

png("Volcano_GSE141910.png", width = 900, height = 700)
plot(
  topGenes$logFC,
  -log10(topGenes$P.Value),
  pch = 20,
  col = ifelse(topGenes$adj.P.Val < 0.05, "red", "grey"),
  main = "Volcano Plot",
  xlab = "log2 Fold Change",
  ylab = "-log10 P-value"
)
abline(h = -log10(0.05), col = "blue", lty = 2)
abline(v = c(-0.1, 0.1), col = "darkgreen", lty = 2)
text(
  top10$logFC,
  -log10(top10$P.Value),
  labels = rownames(top10),
  pos = 3,
  cex = 0.6
)
dev.off()

png("Heatmap_GSE141910.png", width = 900, height = 900)
heatmap(
  dat[top, ],
  Colv = NA,
  scale = "row",
  labCol = FALSE,
  ColSideColors = sample_colors,
  col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  main = "Top Differentially Expressed Genes"
)
dev.off()