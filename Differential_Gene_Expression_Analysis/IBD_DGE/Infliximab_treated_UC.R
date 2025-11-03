#20250922


# DGE for Infliximab treated UC

# Loading packages
library("DESeq2")
library("RColorBrewer")
library("pheatmap")
library("gplots")
library("ggplot2")
library("devtools")

directory <- "/Users/sara/Desktop/Autocore/IBD/counts/"
setwd(directory)

sampleFiles <- list.files(pattern = "_counts.txt$")
sampleNames <- sub("_counts.txt$", "", sampleFiles)

# Keep only important columns
meta <- read.csv("../IBD_Metadata.csv", stringsAsFactors = FALSE)
meta_polished <- meta[, c(1,14,25,31,32)]

# Subset by diagnosis (2nd col) = UC AND treatment (last col) = Infliximab
meta_inflix_UC <- meta_polished[
  toupper(trimws(meta_polished[, 2])) == "UC" &
    toupper(trimws(meta_polished[, ncol(meta_polished)])) == "INFLIXIMAB",
]

# Sample IDs
inflx_UC_ids <- as.character(meta_inflix_UC[, 1])

# Files / names for the UC+Infliximab cohort
Inflix_sampleFiles <- sampleFiles[ sub("_counts.txt$", "", sampleFiles) %in% inflx_UC_ids ]
Inflix_sampleNames <- sub("_counts.txt$", "", Inflix_sampleFiles)

# Conditions (col 3 = Remission / Non-remission) → factor R / NR
status_map <- setNames(as.character(meta_inflix_UC[, 3]), as.character(meta_inflix_UC[, 1]))
Inflix_sampleCondition <- factor(ifelse(status_map[Inflix_sampleNames] == "Remission", "R", "NR"),
                                 levels = c("R","NR"))

## --- DGE for Infliximab + UC cohort 
# Build sample table
sample_table_UC_ifx <- data.frame(
  sampleName = Inflix_sampleNames,
  fileName   = Inflix_sampleFiles,
  condition  = Inflix_sampleCondition,
  stringsAsFactors = FALSE
)

# DESeq2 object and run
ddsHTSeq_UC_ifx <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table_UC_ifx,
                                              design = ~ condition)
colData(ddsHTSeq_UC_ifx)$condition <- factor(colData(ddsHTSeq_UC_ifx)$condition,
                                             levels = c("NR","R"))

dds_UC_ifx <- DESeq(ddsHTSeq_UC_ifx)

# QC: transformed counts, correlation heatmap, PCA
vst_UC_ifx <- varianceStabilizingTransformation(dds_UC_ifx, blind = TRUE)
vst_mat_UC_ifx <- assay(vst_UC_ifx)
vst_cor_UC_ifx <- cor(vst_mat_UC_ifx)

pheatmap(vst_cor_UC_ifx)
plotPCA(vst_UC_ifx, intgroup = "condition")

# rerun the analysis without outlier
if (any(grepl("83257$", Inflix_sampleNames))) {
  keep_ix <- !grepl("83257$", Inflix_sampleNames)
  
  Inflix_sampleNames     <- Inflix_sampleNames[keep_ix]
  Inflix_sampleFiles     <- Inflix_sampleFiles[keep_ix]
  Inflix_sampleCondition <- Inflix_sampleCondition[keep_ix]
  
  sample_table_UC_ifx <- data.frame(
    sampleName = Inflix_sampleNames,
    fileName   = Inflix_sampleFiles,
    condition  = Inflix_sampleCondition,
    stringsAsFactors = FALSE
  )
  
  ddsHTSeq_UC_ifx <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table_UC_ifx,
                                                design = ~ condition)
  colData(ddsHTSeq_UC_ifx)$condition <- factor(colData(ddsHTSeq_UC_ifx)$condition,
                                               levels = c("NR","R"))
  dds_UC_ifx <- DESeq(ddsHTSeq_UC_ifx)
  
  vst_UC_ifx     <- varianceStabilizingTransformation(dds_UC_ifx, blind = TRUE)
  vst_mat_UC_ifx <- assay(vst_UC_ifx)
  vst_cor_UC_ifx <- cor(vst_mat_UC_ifx)
  
  pheatmap(vst_cor_UC_ifx)
  plotPCA(vst_UC_ifx, intgroup = "condition")
}
## PCA seems to be okish, I will proceed with the rest of the analysis.

# Convert Ensembl IDs (human) to HGNC symbols for the UC + Infliximab cohort
library(AnnotationDbi)
library(org.Hs.eg.db)

ensembl_gene_ids_UC_ifx <- sub("\\.\\d+$", "", rownames(dds_UC_ifx))  # remove Ensembl version suffixes
hgnc_symbols_UC_ifx <- mapIds(org.Hs.eg.db,
                              keys      = ensembl_gene_ids_UC_ifx,
                              column    = "SYMBOL",
                              keytype   = "ENSEMBL",
                              multiVals = "first")

mcols(dds_UC_ifx)$hgnc_symbol <- hgnc_symbols_UC_ifx

normalized_counts_UC_ifx <- counts(dds_UC_ifx, normalized = TRUE)
rownames(normalized_counts_UC_ifx) <- ifelse(!is.na(hgnc_symbols_UC_ifx) & hgnc_symbols_UC_ifx != "",
                                             hgnc_symbols_UC_ifx,
                                             rownames(dds_UC_ifx))
View(normalized_counts_UC_ifx)

# Results table (UC + Infliximab) with HGNC symbols for labeling
res_UC_ifx <- results(dds_UC_ifx)
res_UC_ifx$hgnc_symbol <- hgnc_symbols_UC_ifx[rownames(res_UC_ifx)]

res_UC_ifx
summary(res_UC_ifx)

# Filter by adjusted p-value
res05_UC_ifx <- subset(res_UC_ifx, !is.na(padj) & padj < 0.05)
head(res05_UC_ifx)

# Order by padj
resOrdered_UC_ifx <- res05_UC_ifx[order(res05_UC_ifx$padj), ]

# Save (includes symbol column)
write.csv(as.data.frame(resOrdered_UC_ifx), "results-DESeq2-Infliximab_UC.csv", row.names = TRUE)

# MA plot
plotMA(res_UC_ifx, main = "Infliximab_UC_RNA-seq")

# EnhancedVolcano with gene symbols as point labels
library(ggrepel)
library(EnhancedVolcano)

gene_labels_UC_ifx <- ifelse(is.na(res_UC_ifx$hgnc_symbol) | res_UC_ifx$hgnc_symbol == "",
                             rownames(res_UC_ifx), res_UC_ifx$hgnc_symbol)

EnhancedVolcano(as.data.frame(res_UC_ifx),
                lab = gene_labels_UC_ifx,
                x = "log2FoldChange",
                y = "pvalue",
                pCutoff = 1e-3,
                xlim = c(-5, 5))

############
## -------- AutoCore enrichment for Infliximab + UC (SYMBOL-based) --------
library(dplyr)
library(tidyr)
library(clusterProfiler)

# --- Paths ---
autocore_clusters_path <- "/Users/sara/Desktop/Autocore/IBD/Enrichment_Test/Modified_Autocore_Clusters_GeneSymbols copy.txt"
pc_universe_path       <- "/Users/sara/Desktop/Autocore/IBD/Enrichment_Test/symbols2entrezIDs (1).txt"

# 1) DEGs in SYMBOL space (padj < 0.05; adjust if desired)
res_tbl <- as.data.frame(res_UC_ifx) |> tibble::rownames_to_column("ensembl")
stopifnot("hgnc_symbol" %in% names(res_tbl))
res_tbl <- dplyr::rename(res_tbl, symbol = hgnc_symbol)

deg_all_symbols  <- res_tbl |> dplyr::filte_
