# =========================
# FRESH KEGG ANALYSIS
# Female DCM vs Male DCM
# GSE141910
# =========================

library(data.table)
library(GEOquery)
library(Biobase)
library(limma)
library(KEGGREST)
library(ggplot2)

graphics.off()

# =========================
# OUTPUT FOLDER
# =========================

output_folder <- "Fresh_KEGG_Output"
if (!dir.exists(output_folder)) dir.create(output_folder)

# =========================
# LOAD RNA-SEQ COUNT DATA
# =========================

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

cat("Original count matrix dimensions:\n")
print(dim(tbl))

# Filter low-expression genes
keep <- rowSums(tbl >= 10) >= 2
tbl <- tbl[keep, ]

cat("Filtered count matrix dimensions:\n")
print(dim(tbl))

# Log-transform
dat <- log10(tbl + 1)

# =========================
# LOAD METADATA
# =========================

gset <- getGEO("GSE141910", GSEMatrix = TRUE)
geo <- gset[[1]]
pd <- pData(geo)

sample_ids <- colnames(dat)
match_index <- match(sample_ids, pd$geo_accession)
pd <- pd[match_index, ]

# Combine all metadata text so we can detect sex and condition
all_text <- apply(pd, 1, function(x) paste(x, collapse = " | "))
all_text <- tolower(all_text)

# Detect sex
sex <- rep(NA, length(all_text))
sex[grepl("\\bfemale\\b", all_text)] <- "Female"
sex[grepl("\\bmale\\b", all_text) & is.na(sex)] <- "Male"

# Detect disease condition
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

# Keep only samples with known sex and condition
keep_samples <- !is.na(meta$Group)
dat_sex <- dat[, keep_samples]
meta_sex <- meta[keep_samples, ]

cat("\nSample counts by group:\n")
print(table(meta_sex$Group))

# Save metadata
write.csv(
  meta_sex,
  file.path(output_folder, "metadata_used_for_KEGG.csv"),
  row.names = FALSE
)

# =========================
# DIFFERENTIAL EXPRESSION
# Female DCM vs Male DCM
# =========================

dcm_only <- meta_sex$Condition == "DCM"

dat_dcm <- dat_sex[, dcm_only]
meta_dcm <- meta_sex[dcm_only, ]

meta_dcm$Sex <- factor(meta_dcm$Sex, levels = c("Male", "Female"))

cat("\nDCM sample counts by sex:\n")
print(table(meta_dcm$Sex))

design <- model.matrix(~ meta_dcm$Sex)

fit <- lmFit(dat_dcm, design)
fit <- eBayes(fit)

topGenes <- topTable(fit, coef = 2, number = Inf)

write.csv(
  topGenes,
  file.path(output_folder, "DE_Female_DCM_vs_Male_DCM.csv"),
  row.names = TRUE
)

cat("\nDifferential expression completed.\n")
cat("Total genes tested:", nrow(topGenes), "\n")

# =========================
# SELECT GENES FOR KEGG
# =========================

# Main cutoff
sig_genes <- topGenes[topGenes$adj.P.Val < 0.05, ]

# If too few genes pass adjusted p-value, use nominal p-value for exploratory KEGG
if (nrow(sig_genes) < 10) {
  message("Fewer than 10 genes passed adj.P.Val < 0.05. Using P.Value < 0.05 for exploratory KEGG.")
  sig_genes <- topGenes[topGenes$P.Value < 0.05, ]
  gene_cutoff_used <- "P.Value < 0.05"
} else {
  gene_cutoff_used <- "adj.P.Val < 0.05"
}

sig_entrez <- unique(as.character(rownames(sig_genes)))
universe_entrez <- unique(as.character(rownames(topGenes)))

# Keep only numeric Entrez-style IDs
sig_entrez <- sig_entrez[grepl("^[0-9]+$", sig_entrez)]
universe_entrez <- universe_entrez[grepl("^[0-9]+$", universe_entrez)]

cat("\nGenes used for KEGG:", length(sig_entrez), "\n")
cat("Universe genes:", length(universe_entrez), "\n")
cat("Cutoff used:", gene_cutoff_used, "\n")

write.csv(
  data.frame(Entrez_ID = sig_entrez),
  file.path(output_folder, "genes_used_for_KEGG.csv"),
  row.names = FALSE
)

# =========================
# GET HUMAN KEGG PATHWAYS
# =========================

cat("\nDownloading KEGG pathway annotations...\n")

kegg_pathways <- keggList("pathway", "hsa")
kegg_links <- keggLink("pathway", "hsa")

# Convert KEGG IDs like hsa:1234 to 1234
kegg_gene_ids <- sub("hsa:", "", names(kegg_links))

pathway_ids <- sub("path:", "", as.character(kegg_links))

kegg_map <- data.frame(
  GeneID = kegg_gene_ids,
  PathwayID = pathway_ids,
  stringsAsFactors = FALSE
)

# Add pathway names
pathway_names <- data.frame(
  PathwayID = sub("path:", "", names(kegg_pathways)),
  PathwayName = as.character(kegg_pathways),
  stringsAsFactors = FALSE
)

pathway_names$PathwayName <- sub(" - Homo sapiens \\(human\\)", "", pathway_names$PathwayName)

# Restrict universe to genes that are in KEGG
kegg_universe <- intersect(universe_entrez, unique(kegg_map$GeneID))
sig_in_kegg <- intersect(sig_entrez, kegg_universe)

cat("Universe genes found in KEGG:", length(kegg_universe), "\n")
cat("Significant genes found in KEGG:", length(sig_in_kegg), "\n")

# =========================
# HYPERGEOMETRIC KEGG ENRICHMENT
# =========================

results <- data.frame()

for (pid in unique(kegg_map$PathwayID)) {
  
  pathway_genes <- unique(kegg_map$GeneID[kegg_map$PathwayID == pid])
  pathway_genes <- intersect(pathway_genes, kegg_universe)
  
  overlap_genes <- intersect(sig_in_kegg, pathway_genes)
  
  k <- length(overlap_genes)                 # significant genes in pathway
  M <- length(pathway_genes)                 # genes in pathway
  N <- length(kegg_universe)                 # all genes in KEGG universe
  n <- length(sig_in_kegg)                   # significant genes in KEGG
  
  if (M > 0 && n > 0 && k > 0) {
    
    pval <- phyper(
      q = k - 1,
      m = M,
      n = N - M,
      k = n,
      lower.tail = FALSE
    )
    
    results <- rbind(
      results,
      data.frame(
        PathwayID = pid,
        PathwaySize = M,
        SignificantGenesInPathway = k,
        SignificantGenesTotal = n,
        GeneRatio = k / n,
        BackgroundRatio = M / N,
        P.Value = pval,
        Genes = paste(overlap_genes, collapse = "/"),
        stringsAsFactors = FALSE
      )
    )
  }
}

if (nrow(results) > 0) {
  
  results$adj.P.Val <- p.adjust(results$P.Value, method = "BH")
  
  results <- merge(results, pathway_names, by = "PathwayID", all.x = TRUE)
  results <- results[order(results$P.Value), ]
  
} else {
  
  results <- data.frame(
    Message = "No KEGG pathways had overlap with the significant gene list."
  )
}

write.csv(
  results,
  file.path(output_folder, "Fresh_KEGG_Enrichment_Results.csv"),
  row.names = FALSE
)

# =========================
# MAKE KEGG DOTPLOT
# =========================

if ("P.Value" %in% colnames(results) && nrow(results) > 0) {
  
  plot_df <- results[order(results$P.Value), ]
  plot_df <- head(plot_df, 15)
  
  plot_df$PathwayName <- factor(
    plot_df$PathwayName,
    levels = rev(plot_df$PathwayName)
  )
  
  plot_title <- ifelse(
    any(plot_df$adj.P.Val < 0.05),
    "KEGG Pathway Enrichment: Female DCM vs Male DCM",
    "Top KEGG Pathways: Female DCM vs Male DCM"
  )
  
  plot_subtitle <- ifelse(
    any(plot_df$adj.P.Val < 0.05),
    paste0("Gene cutoff: ", gene_cutoff_used, " | FDR-significant pathways detected"),
    paste0("Gene cutoff: ", gene_cutoff_used, " | No pathways pass FDR < 0.05")
  )
  
  kegg_plot <- ggplot(
    plot_df,
    aes(
      x = GeneRatio,
      y = PathwayName,
      size = SignificantGenesInPathway,
      color = adj.P.Val
    )
  ) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(
      low = "red",
      high = "blue",
      name = "Adjusted\nP-value"
    ) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Gene Ratio",
      y = "KEGG Pathway",
      size = "Gene Count"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      axis.text.y = element_text(size = 10),
      legend.position = "right"
    )
  
  print(kegg_plot)
  
  ggsave(
    filename = file.path(output_folder, "Fresh_KEGG_Dotplot_Female_DCM_vs_Male_DCM.png"),
    plot = kegg_plot,
    width = 11,
    height = 7,
    dpi = 300
  )
  
  # Also make a clean barplot
  bar_plot <- ggplot(
    plot_df,
    aes(
      x = SignificantGenesInPathway,
      y = PathwayName,
      fill = adj.P.Val
    )
  ) +
    geom_col() +
    scale_fill_gradient(
      low = "red",
      high = "blue",
      name = "Adjusted\nP-value"
    ) +
    labs(
      title = "KEGG Pathway Gene Counts",
      subtitle = plot_subtitle,
      x = "Number of significant genes in pathway",
      y = "KEGG Pathway"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      axis.text.y = element_text(size = 10),
      legend.position = "right"
    )
  
  print(bar_plot)
  
  ggsave(
    filename = file.path(output_folder, "Fresh_KEGG_Barplot_Female_DCM_vs_Male_DCM.png"),
    plot = bar_plot,
    width = 11,
    height = 7,
    dpi = 300
  )
  
} else {
  
  message("No KEGG plot created because no pathways overlapped with the gene list.")
}

# =========================
# WRITE SUMMARY FILE
# =========================

summary_lines <- c(
  "Fresh KEGG Analysis Summary",
  "Comparison: Female DCM vs Male DCM",
  paste("Gene cutoff used:", gene_cutoff_used),
  paste("Total genes tested:", nrow(topGenes)),
  paste("Genes selected for KEGG:", length(sig_entrez)),
  paste("Universe genes:", length(universe_entrez)),
  paste("Universe genes found in KEGG:", length(kegg_universe)),
  paste("Significant genes found in KEGG:", length(sig_in_kegg))
)

if ("adj.P.Val" %in% colnames(results)) {
  summary_lines <- c(
    summary_lines,
    paste("KEGG pathways with overlap:", nrow(results)),
    paste("FDR-significant KEGG pathways:", sum(results$adj.P.Val < 0.05, na.rm = TRUE)),
    paste("Nominally significant KEGG pathways:", sum(results$P.Value < 0.05, na.rm = TRUE))
  )
} else {
  summary_lines <- c(
    summary_lines,
    "KEGG pathways with overlap: 0"
  )
}

writeLines(
  summary_lines,
  file.path(output_folder, "Fresh_KEGG_Summary.txt")
)

cat("\nDONE. Check this folder:\n")
cat(output_folder, "\n")
cat("\nFiles created:\n")
print(list.files(output_folder))