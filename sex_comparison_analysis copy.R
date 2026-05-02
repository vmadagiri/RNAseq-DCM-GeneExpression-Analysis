# -------------------------
# Libraries
# -------------------------

library(data.table)
library(GEOquery)
library(Biobase)
library(limma)
library(umap)

graphics.off()

# -------------------------
# Load data
# -------------------------

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

keep <- rowSums(tbl >= 10) >= 2
tbl <- tbl[keep, ]

dat <- log10(tbl + 1)

# -------------------------
# Metadata
# -------------------------

gset <- getGEO("GSE141910", GSEMatrix = TRUE)
geo <- gset[[1]]
pd <- pData(geo)

sample_ids <- colnames(dat)

match_index <- match(sample_ids, pd$geo_accession)
pd <- pd[match_index, ]

all_text <- apply(pd, 1, function(x) paste(x, collapse = " | "))
all_text <- tolower(all_text)

sex <- rep(NA, length(all_text))
sex[grepl("\\bfemale\\b", all_text)] <- "Female"
sex[grepl("\\bmale\\b", all_text) & is.na(sex)] <- "Male"

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

meta_sex$Sex <- factor(meta_sex$Sex, levels = c("Male", "Female"))
meta_sex$Condition <- factor(meta_sex$Condition, levels = c("Control", "DCM"))

# colors that are easier to read
group_colors <- c(
  "Male_Control" = "blue",
  "Female_Control" = "deeppink",
  "Male_DCM" = "darkgreen",
  "Female_DCM" = "gold2"
)

sex_colors <- c(
  "Male" = "blue",
  "Female" = "deeppink"
)

output_folder <- "Sex_Comparison_Clean_Final"

if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# -------------------------
# PCA
# -------------------------

pca <- prcomp(t(dat_sex), scale. = TRUE)

pc1_var <- round(100 * summary(pca)$importance[2, 1], 1)
pc2_var <- round(100 * summary(pca)$importance[2, 2], 1)

png(
  file.path(output_folder, "01_PCA_Sex_Comparison.png"),
  width = 1800,
  height = 1100,
  res = 180
)

par(mar = c(5, 5, 4, 11))

plot(
  pca$x[, 1],
  pca$x[, 2],
  col = group_colors[as.character(meta_sex$Group)],
  pch = 19,
  cex = 0.85,
  xlab = paste0("PC1 (", pc1_var, "%)"),
  ylab = paste0("PC2 (", pc2_var, "%)"),
  main = "PCA Plot - Sex and Disease Groups"
)

par(xpd = TRUE)

legend(
  "topright",
  inset = c(-0.30, 0),
  legend = levels(meta_sex$Group),
  col = group_colors[levels(meta_sex$Group)],
  pch = 19,
  cex = 0.85,
  bty = "n"
)

par(xpd = FALSE)

dev.off()

# -------------------------
# UMAP
# -------------------------

dat_umap <- dat_sex[!duplicated(dat_sex), ]

set.seed(123)

ump <- umap(
  t(dat_umap),
  n_neighbors = 20,
  min_dist = 0.2,
  random_state = 123
)

png(
  file.path(output_folder, "02_UMAP_Sex_Comparison.png"),
  width = 1800,
  height = 1100,
  res = 180
)

par(mar = c(5, 5, 4, 11))

plot(
  ump$layout[, 1],
  ump$layout[, 2],
  col = group_colors[as.character(meta_sex$Group)],
  pch = 19,
  cex = 0.9,
  xlab = "UMAP1",
  ylab = "UMAP2",
  main = "UMAP Plot - Sex and Disease Groups"
)

par(xpd = TRUE)

legend(
  "topright",
  inset = c(-0.30, 0),
  legend = levels(meta_sex$Group),
  col = group_colors[levels(meta_sex$Group)],
  pch = 19,
  cex = 0.85,
  bty = "n"
)

par(xpd = FALSE)

dev.off()

# -------------------------
# PC1 boxplot
# -------------------------

pc1_df <- data.frame(
  PC1 = pca$x[, 1],
  Group = meta_sex$Group
)

png(
  file.path(output_folder, "03_PC1_Boxplot.png"),
  width = 1800,
  height = 1100,
  res = 180
)

par(mar = c(8, 5, 4, 11))

boxplot(
  PC1 ~ Group,
  data = pc1_df,
  col = group_colors[levels(meta_sex$Group)],
  border = "black",
  las = 2,
  outline = FALSE,
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
  cex = 0.45,
  col = rgb(0, 0, 0, 0.5),
  add = TRUE
)

par(xpd = TRUE)

legend(
  "topright",
  inset = c(-0.30, 0),
  legend = levels(meta_sex$Group),
  fill = group_colors[levels(meta_sex$Group)],
  cex = 0.85,
  bty = "n"
)

par(xpd = FALSE)

dev.off()

# -------------------------
# Volcano function
# -------------------------

make_volcano <- function(dat_input, meta_input, title_name, file_name) {
  
  meta_input$Sex <- factor(meta_input$Sex, levels = c("Male", "Female"))
  
  design <- model.matrix(~ meta_input$Sex)
  
  fit <- lmFit(dat_input, design)
  fit <- eBayes(fit)
  
  topGenes <- topTable(fit, coef = 2, number = Inf)
  
  volcano_df <- topGenes
  
  volcano_df$P.Value[is.na(volcano_df$P.Value)] <- 1
  volcano_df$P.Value[volcano_df$P.Value <= 0] <- 1e-300
  
  volcano_df$negLogP <- -log10(volcano_df$P.Value)
  
  # this keeps the graph from getting stretched by extreme genes
  volcano_df$plotNegLogP <- volcano_df$negLogP
  volcano_df$plotNegLogP[volcano_df$plotNegLogP > 25] <- 25
  
  p_cutoff <- 0.05
  logfc_cutoff <- 0.20
  
  volcano_df$Significance <- "Not Significant"
  
  volcano_df$Significance[
    volcano_df$P.Value < p_cutoff &
      volcano_df$logFC > logfc_cutoff
  ] <- "Female Higher"
  
  volcano_df$Significance[
    volcano_df$P.Value < p_cutoff &
      volcano_df$logFC < -logfc_cutoff
  ] <- "Male Higher"
  
  volcano_df$Color <- "grey80"
  volcano_df$Color[volcano_df$Significance == "Female Higher"] <- "deeppink"
  volcano_df$Color[volcano_df$Significance == "Male Higher"] <- "blue"
  
  png(
    file.path(output_folder, file_name),
    width = 1800,
    height = 1100,
    res = 180
  )
  
  par(mar = c(5, 5, 4, 11))
  
  plot(
    volcano_df$logFC,
    volcano_df$plotNegLogP,
    col = adjustcolor(volcano_df$Color, alpha.f = 0.55),
    pch = 19,
    cex = 0.55,
    xlim = c(-1.5, 1.5),
    ylim = c(0, 26),
    xlab = "log Fold Change (Female vs Male)",
    ylab = "-log10 P-value",
    main = title_name
  )
  
  abline(h = -log10(p_cutoff), col = "darkgreen", lty = 2, lwd = 1.5)
  abline(v = c(-logfc_cutoff, logfc_cutoff), col = "darkgreen", lty = 2, lwd = 1.5)
  
  par(xpd = TRUE)
  
  legend(
    "topright",
    inset = c(-0.30, 0),
    legend = c("Female Higher", "Male Higher", "Not Significant"),
    col = c("deeppink", "blue", "grey80"),
    pch = 19,
    cex = 0.85,
    bty = "n"
  )
  
  par(xpd = FALSE)
  
  dev.off()
  
  print(title_name)
  print(table(volcano_df$Significance))
  
  return(volcano_df)
}

# -------------------------
# DCM volcano
# -------------------------

dcm_only <- meta_sex$Condition == "DCM"

dat_dcm <- dat_sex[, dcm_only]
meta_dcm <- meta_sex[dcm_only, ]

volcano_dcm_df <- make_volcano(
  dat_dcm,
  meta_dcm,
  "Volcano Plot - Female vs Male DCM",
  "04_Volcano_DCM_Clean.png"
)

# -------------------------
# Control volcano
# -------------------------

control_only <- meta_sex$Condition == "Control"

dat_control <- dat_sex[, control_only]
meta_control <- meta_sex[control_only, ]

volcano_control_df <- make_volcano(
  dat_control,
  meta_control,
  "Volcano Plot - Female vs Male Control",
  "05_Volcano_Control_Clean.png"
)

# -------------------------
# Heatmap using top sex-difference genes
# -------------------------

# use strongest DCM sex-difference genes for heatmap
heat_genes <- rownames(volcano_dcm_df[order(volcano_dcm_df$P.Value), ])[1:50]

heat_data <- dat_dcm[heat_genes, ]

heat_data_scaled <- t(scale(t(heat_data)))

heat_data_scaled[is.na(heat_data_scaled)] <- 0
heat_data_scaled[heat_data_scaled > 2] <- 2
heat_data_scaled[heat_data_scaled < -2] <- -2

heat_colors <- sex_colors[as.character(meta_dcm$Sex)]

png(
  file.path(output_folder, "06_Heatmap_Top_DCM_Sex_Difference_Genes.png"),
  width = 1800,
  height = 1200,
  res = 180
)

heatmap(
  heat_data_scaled,
  ColSideColors = heat_colors,
  scale = "none",
  labRow = FALSE,
  labCol = FALSE,
  margins = c(5, 5),
  col = colorRampPalette(c("blue", "white", "deeppink"))(75),
  main = "Heatmap - Top DCM Sex Difference Genes"
)

legend(
  "topright",
  legend = c("Male DCM", "Female DCM"),
  fill = c("blue", "deeppink"),
  cex = 0.85,
  bty = "n"
)

dev.off()

# -------------------------
# Save tables
# -------------------------

write.csv(
  volcano_dcm_df,
  file.path(output_folder, "Volcano_DCM_Table.csv")
)

write.csv(
  volcano_control_df,
  file.path(output_folder, "Volcano_Control_Table.csv")
)

write.csv(
  meta_sex,
  file.path(output_folder, "Metadata_Used.csv"),
  row.names = FALSE
)

# -------------------------
# Done
# -------------------------

list.files(output_folder)

print("DONE. Clean plots saved in Sex_Comparison_Clean_Final")

