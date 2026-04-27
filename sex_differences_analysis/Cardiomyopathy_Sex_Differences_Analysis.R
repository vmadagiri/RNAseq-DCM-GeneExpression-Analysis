library(data.table)
library(GEOquery)
library(Biobase)
library(limma)
library(umap)

graphics.off()

# make sure plots go to RStudio Plots pane
if ("RStudioGD" %in% names(grDevices::deviceIsInteractive)) {
  options(device = "RStudioGD")
}

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

group_colors <- c(
  "Male_Control" = "blue",
  "Female_Control" = "deeppink",
  "Male_DCM" = "red",
  "Female_DCM" = "purple"
)

output_folder <- "Sex_Comparison_Plots"
if (!dir.exists(output_folder)) dir.create(output_folder)

# -------------------------
# PCA
# -------------------------

pca <- prcomp(t(dat_sex), scale. = TRUE)
pc1_var <- round(100 * summary(pca)$importance[2, 1], 1)
pc2_var <- round(100 * summary(pca)$importance[2, 2], 1)

par(mar = c(5, 5, 4, 11))
plot(
  pca$x[, 1],
  pca$x[, 2],
  col = group_colors[as.character(meta_sex$Group)],
  pch = 19,
  cex = 0.9,
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
dev.copy(png, file.path(output_folder, "PCA_Sex_Comparison_GSE141910.png"), width = 1800, height = 1000, res = 180)
dev.off()

# -------------------------
# UMAP
# -------------------------

dat2 <- dat_sex[!duplicated(dat_sex), ]
ump <- umap(t(dat2), n_neighbors = 15, random_state = 123)

par(mar = c(5, 5, 4, 11))
plot(
  ump$layout[, 1],
  ump$layout[, 2],
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
dev.copy(png, file.path(output_folder, "UMAP_Sex_Comparison_GSE141910.png"), width = 1800, height = 1000, res = 180)
dev.off()

# -------------------------
# Boxplot
# -------------------------

pc1_df <- data.frame(
  PC1 = pca$x[, 1],
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
dev.copy(png, file.path(output_folder, "PC1_Boxplot_Sex_Comparison_GSE141910.png"), width = 1800, height = 1000, res = 180)
dev.off()

# -------------------------
# Volcano
# -------------------------

dcm_only <- meta_sex$Condition == "DCM"
dat_dcm <- dat_sex[, dcm_only]
meta_dcm <- meta_sex[dcm_only, ]
meta_dcm$Sex <- factor(meta_dcm$Sex, levels = c("Male", "Female"))

design_sex <- model.matrix(~ meta_dcm$Sex)
fit_sex <- lmFit(dat_dcm, design_sex)
fit_sex <- eBayes(fit_sex)
topGenes_sex <- topTable(fit_sex, coef = 2, number = Inf)

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

top_up <- head(volcano_df[volcano_df$Significance == "Upregulated", ][order(volcano_df[volcano_df$Significance == "Upregulated", ]$adj.P.Val), ], 3)
top_down <- head(volcano_df[volcano_df$Significance == "Downregulated", ][order(volcano_df[volcano_df$Significance == "Downregulated", ]$adj.P.Val), ], 3)
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
dev.copy(png, file.path(output_folder, "Volcano_Male_vs_Female_DCM_GSE141910.png"), width = 1800, height = 1000, res = 180)
dev.off()

dev.copy(png, file.path(output_folder, "Volcano_Male_vs_Female_DCM_GSE141910.png"), width = 1800, height = 1000, res = 180)
dev.off()

# -------------------------
# NEW CONTROL VOLCANO
# -------------------------

control_only <- meta_sex$Condition == "Control"
dat_control <- dat_sex[, control_only]
meta_control <- meta_sex[control_only, ]

meta_control$Sex <- factor(meta_control$Sex, levels = c("Male", "Female"))

design_control <- model.matrix(~ meta_control$Sex)
fit_control <- lmFit(dat_control, design_control)
fit_control <- eBayes(fit_control)

topGenes_control <- topTable(fit_control, coef = 2, number = Inf)

volcano_control_df <- topGenes_control
volcano_control_df$adj.P.Val[is.na(volcano_control_df$adj.P.Val)] <- 1
volcano_control_df$adj.P.Val[volcano_control_df$adj.P.Val <= 0] <- 1e-300
volcano_control_df$negLogAdjP <- -log10(volcano_control_df$adj.P.Val)

volcano_control_df$Significance <- "Not Significant"
volcano_control_df$Significance[volcano_control_df$adj.P.Val < 0.05 & volcano_control_df$logFC > 1] <- "Upregulated"
volcano_control_df$Significance[volcano_control_df$adj.P.Val < 0.05 & volcano_control_df$logFC < -1] <- "Downregulated"

volcano_control_df$Color <- "grey"
volcano_control_df$Color[volcano_control_df$Significance == "Upregulated"] <- "red"
volcano_control_df$Color[volcano_control_df$Significance == "Downregulated"] <- "blue"

y_max_control <- max(volcano_control_df$negLogAdjP, na.rm = TRUE)
y_pad_control <- 0.08 * y_max_control

par(mar = c(5, 5, 4, 11))
plot(
  volcano_control_df$logFC,
  volcano_control_df$negLogAdjP,
  col = volcano_control_df$Color,
  pch = 19,
  cex = 0.65,
  ylim = c(0, y_max_control + y_pad_control),
  xlab = "log Fold Change",
  ylab = "-log10 Adjusted P-value",
  main = "Volcano Plot - Female vs Male Control"
)

abline(h = -log10(0.05), col = "darkgreen", lty = 2, lwd = 1.5)
abline(v = c(-1, 1), col = "darkgreen", lty = 2, lwd = 1.5)

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
  file.path(output_folder, "Volcano_Control_Sex_Comparison_GSE141910.png"),
  width = 1800,
  height = 1000,
  res = 180
)
dev.off()

list.files(output_folder)

list.files(output_folder)