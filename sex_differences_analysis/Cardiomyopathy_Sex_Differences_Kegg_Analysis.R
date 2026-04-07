library(data.table)
library(GEOquery)
library(Biobase)
library(limma)
library(umap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

graphics.off()

# get the data
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"

path <- paste(
  urld,
  "acc=GSE141910",
  "file=GSE141910_raw_counts_GRCh38.p13_NCBI.tsv.gz",
  sep = "&"
)

tbl <- as.matrix(
  fread(path, header = TRUE, colClasses = "integer"),
  rownames = 1
)

# keep genes with at least a little expression
keep <- rowSums(tbl >= 10) >= 2
tbl <- tbl[keep, ]

# log transform
dat <- log10(tbl + 1)

# get metadata
gset <- getGEO("GSE141910", GSEMatrix = TRUE)
geo <- gset[[1]]
pd <- pData(geo)

sample_ids <- colnames(dat)
match_index <- match(sample_ids, pd$geo_accession)
pd <- pd[match_index, ]

all_text <- apply(pd, 1, function(x) paste(x, collapse = " | "))
all_text <- tolower(all_text)

# figure out sex
sex <- rep(NA, length(all_text))
sex[grepl("\\bfemale\\b", all_text)] <- "Female"
sex[grepl("\\bmale\\b", all_text) & is.na(sex)] <- "Male"

# figure out condition
condition <- rep(NA, length(all_text))
condition[grepl("\\bdcm\\b", all_text)] <- "DCM"
condition[
  grepl("non-failing", all_text) |
    grepl("non failing", all_text) |
    grepl("\\bcontrol\\b", all_text) |
    grepl("\\bnf\\b", all_text)
] <- "Control"

meta <- data.frame(
  Sample = sample_ids,
  Sex = sex,
  Condition = condition,
  stringsAsFactors = FALSE
)

meta$Group <- paste(meta$Sex, meta$Condition, sep = "_")
meta$Group[is.na(meta$Sex) | is.na(meta$Condition)] <- NA

keep_samples <- !is.na(meta$Group)
dat_sex <- dat[, keep_samples]
meta_sex <- meta[keep_samples, ]

meta_sex$Group <- factor(
  meta_sex$Group,
  levels = c("Male_Control", "Female_Control", "Male_DCM", "Female_DCM")
)

group_colors <- c(
  "Male_Control" = "blue",
  "Female_Control" = "deeppink",
  "Male_DCM" = "red",
  "Female_DCM" = "purple"
)

output_folder <- "Sex_Comparison_Plots"
if (!dir.exists(output_folder)) dir.create(output_folder)

# PCA
pca <- prcomp(t(dat_sex), scale. = TRUE)

pc1_var <- round(100 * summary(pca)$importance[2, 1], 1)
pc2_var <- round(100 * summary(pca)$importance[2, 2], 1)

par(mar = c(5, 5, 4, 11))
plot(
  pca$x[,1],
  pca$x[,2],
  col = group_colors[as.character(meta_sex$Group)],
  pch = 19,
  cex = 1,
  xlab = paste0("PC1 (", pc1_var, "%)"),
  ylab = paste0("PC2 (", pc2_var, "%)"),
  main = "PCA Plot - Sex Comparison"
)

par(xpd = TRUE)
legend(
  "topright",
  inset = c(-0.28, 0),
  legend = levels(meta_sex$Group),
  col = group_colors[levels(meta_sex$Group)],
  pch = 19,
  cex = 0.85,
  bty = "n"
)
par(xpd = FALSE)

dev.copy(
  png,
  file.path(output_folder, "PCA_Sex_Comparison_GSE141910.png"),
  width = 1800,
  height = 1000,
  res = 180
)
dev.off()

# UMAP
dat2 <- dat_sex[!duplicated(dat_sex), ]
ump <- umap(t(dat2), n_neighbors = 15, random_state = 123)

par(mar = c(5, 5, 4, 11))
plot(
  ump$layout[,1],
  ump$layout[,2],
  col = group_colors[as.character(meta_sex$Group)],
  pch = 19,
  cex = 0.9,
  xlab = "UMAP1",
  ylab = "UMAP2",
  main = "UMAP Plot - Sex Comparison"
)

par(xpd = TRUE)
legend(
  "topright",
  inset = c(-0.28, 0),
  legend = levels(meta_sex$Group),
  col = group_colors[levels(meta_sex$Group)],
  pch = 19,
  cex = 0.85,
  bty = "n"
)
par(xpd = FALSE)

dev.copy(
  png,
  file.path(output_folder, "UMAP_Sex_Comparison_GSE141910.png"),
  width = 1800,
  height = 1000,
  res = 180
)
dev.off()

# boxplot for PC1
pc1_df <- data.frame(
  PC1 = pca$x[,1],
  Group = meta_sex$Group
)

y_range <- range(pc1_df$PC1, na.rm = TRUE)
y_pad <- 0.08 * diff(y_range)

par(mar = c(8, 5, 4, 11))
boxplot(
  PC1 ~ Group,
  data = pc1_df,
  col = group_colors[levels(meta_sex$Group)],
  border = "black",
  las = 2,
  outline = FALSE,
  ylim = c(y_range[1] - y_pad, y_range[2] + y_pad),
  ylab = "PC1 Score",
  xlab = "",
  main = "PC1 Distribution by Group"
)

stripchart(
  PC1 ~ Group,
  data = pc1_df,
  vertical = TRUE,
  method = "jitter",
  jitter = 0.12,
  pch = 19,
  cex = 0.5,
  col = rgb(0, 0, 0, 0.6),
  add = TRUE
)

par(xpd = TRUE)
legend(
  "topright",
  inset = c(-0.28, 0),
  legend = levels(meta_sex$Group),
  fill = group_colors[levels(meta_sex$Group)],
  cex = 0.85,
  bty = "n"
)
par(xpd = FALSE)

dev.copy(
  png,
  file.path(output_folder, "PC1_Boxplot_Sex_Comparison_GSE141910.png"),
  width = 1800,
  height = 1000,
  res = 180
)
dev.off()

# only looking at DCM for the sex comparison volcano
dcm_only <- meta_sex$Condition == "DCM"
dat_dcm <- dat_sex[, dcm_only]
meta_dcm <- meta_sex[dcm_only, ]

meta_dcm$Sex <- factor(meta_dcm$Sex, levels = c("Male", "Female"))

design_sex <- model.matrix(~ meta_dcm$Sex)
fit_sex <- lmFit(dat_dcm, design_sex)
fit_sex <- eBayes(fit_sex)

topGenes_sex <- topTable(fit_sex, coef = 2, number = Inf)

# KEGG
sig_genes_kegg <- topGenes_sex[topGenes_sex$adj.P.Val < 0.05, ]

cat("Significant genes used for KEGG:", nrow(sig_genes_kegg), "\n")

kegg_gene_ids <- rownames(sig_genes_kegg)
kegg_gene_ids <- as.character(kegg_gene_ids)
kegg_gene_ids <- trimws(kegg_gene_ids)
kegg_gene_ids <- kegg_gene_ids[kegg_gene_ids != ""]
kegg_gene_ids <- unique(kegg_gene_ids)

all_gene_ids <- rownames(topGenes_sex)
all_gene_ids <- as.character(all_gene_ids)
all_gene_ids <- trimws(all_gene_ids)
all_gene_ids <- all_gene_ids[all_gene_ids != ""]
all_gene_ids <- unique(all_gene_ids)

cat("Gene IDs going into KEGG:", length(kegg_gene_ids), "\n")

kegg_df <- data.frame()

if (length(kegg_gene_ids) > 0) {
  kegg_results <- tryCatch(
    enrichKEGG(
      gene = kegg_gene_ids,
      organism = "hsa",
      universe = all_gene_ids,
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.2
    ),
    error = function(e) NULL
  )
  
  if (!is.null(kegg_results)) {
    kegg_df <- as.data.frame(kegg_results)
  }
}

cat("Enriched KEGG pathways found:", nrow(kegg_df), "\n")

# just made this so I had something visual for KEGG
counts <- c(nrow(sig_genes_kegg), max(1, nrow(kegg_df)))

png(
  file.path(output_folder, "KEGG_Summary_Barplot.png"),
  width = 1200,
  height = 800
)

barplot(
  counts,
  names.arg = c("Significant Genes", "Enriched Pathways"),
  main = "KEGG Analysis Summary",
  ylab = "Count",
  ylim = c(0, max(counts) * 1.2),
  col = c("green", "orange")
)

text(
  x = c(0.7, 1.9),
  y = counts,
  labels = c(nrow(sig_genes_kegg), nrow(kegg_df)),
  pos = 3,
  cex = 1
)

dev.off()

if (nrow(kegg_df) > 0) {
  write.csv(
    kegg_df,
    file.path(output_folder, "KEGG_Sex_Comparison_Results.csv"),
    row.names = FALSE
  )
  
  png(
    file.path(output_folder, "KEGG_Sex_Comparison_Dotplot.png"),
    width = 1800,
    height = 1000,
    res = 180
  )
  print(dotplot(kegg_results, showCategory = 10))
  dev.off()
  
  png(
    file.path(output_folder, "KEGG_Sex_Comparison_Barplot.png"),
    width = 1800,
    height = 1000,
    res = 180
  )
  print(barplot(kegg_results, showCategory = 10, drop = TRUE))
  dev.off()
  
  writeLines(
    c(
      "KEGG pathway analysis ran successfully.",
      paste("Significant genes tested:", nrow(sig_genes_kegg)),
      paste("Gene IDs used:", length(kegg_gene_ids)),
      paste("Enriched pathways found:", nrow(kegg_df))
    ),
    file.path(output_folder, "KEGG_Sex_Comparison_Summary.txt")
  )
  
  message("KEGG finished and the files were saved.")
  
} else {
  write.csv(
    data.frame(
      Message = "KEGG was run but no enriched pathways were found.",
      Significant_Genes_Tested = nrow(sig_genes_kegg),
      Gene_IDs_Used = length(kegg_gene_ids)
    ),
    file.path(output_folder, "KEGG_Sex_Comparison_Results.csv"),
    row.names = FALSE
  )
  
  writeLines(
    c(
      "KEGG pathway analysis was done.",
      paste("Significant genes tested:", nrow(sig_genes_kegg)),
      paste("Gene IDs used:", length(kegg_gene_ids)),
      "No enriched KEGG pathways were found."
    ),
    file.path(output_folder, "KEGG_Sex_Comparison_Summary.txt")
  )
  
  message("KEGG ran, but there were no enriched pathways.")
}

# volcano plot
volcano_df <- topGenes_sex

volcano_df$adj.P.Val[is.na(volcano_df$adj.P.Val)] <- 1
volcano_df$adj.P.Val[volcano_df$adj.P.Val <= 0] <- 1e-300
volcano_df$negLogAdjP <- -log10(volcano_df$adj.P.Val)

volcano_df$Significance <- "Not Significant"
volcano_df$Significance[volcano_df$adj.P.Val < 0.05 & volcano_df$logFC > 1] <- "Upregulated"
volcano_df$Significance[volcano_df$adj.P.Val < 0.05 & volcano_df$logFC < -1] <- "Downregulated"

volcano_df$Color <- "grey"
volcano_df$Color[volcano_df$Significance == "Upregulated"] <- "red"
volcano_df$Color[volcano_df$Significance == "Downregulated"] <- "blue"

top_up <- head(
  volcano_df[volcano_df$Significance == "Upregulated", ][
    order(volcano_df[volcano_df$Significance == "Upregulated", ]$adj.P.Val),
  ],
  3
)

top_down <- head(
  volcano_df[volcano_df$Significance == "Downregulated", ][
    order(volcano_df[volcano_df$Significance == "Downregulated", ]$adj.P.Val),
  ],
  3
)

top_label_df <- rbind(top_up, top_down)

y_max <- max(volcano_df$negLogAdjP, na.rm = TRUE)
y_pad <- 0.08 * y_max

par(mar = c(5, 5, 4, 11))
plot(
  volcano_df$logFC,
  volcano_df$negLogAdjP,
  col = volcano_df$Color,
  pch = 19,
  cex = 0.65,
  ylim = c(0, y_max + y_pad),
  xlab = "log Fold Change",
  ylab = "-log10 Adjusted P-value",
  main = "Volcano Plot - Female vs Male DCM"
)

abline(h = -log10(0.05), col = "darkgreen", lty = 2, lwd = 1.5)
abline(v = c(-1, 1), col = "darkgreen", lty = 2, lwd = 1.5)

if (nrow(top_label_df) > 0) {
  text(
    top_label_df$logFC,
    top_label_df$negLogAdjP,
    labels = rownames(top_label_df),
    pos = 3,
    cex = 0.5,
    offset = 0.3
  )
}

par(xpd = TRUE)
legend(
  "topright",
  inset = c(-0.28, 0),
  legend = c("Upregulated", "Downregulated", "Not Significant"),
  col = c("red", "blue", "grey"),
  pch = 19,
  cex = 0.85,
  bty = "n"
)
par(xpd = FALSE)

dev.copy(
  png,
  file.path(output_folder, "Volcano_Sex_Comparison_GSE141910.png"),
  width = 1800,
  height = 1000,
  res = 180
)
dev.off()