
#20251001
# This R snippet is modified by ChatGPT
# Prompt: DGE for Infliximab treated CD

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

# Subset by diagnosis (2nd col) = CD AND treatment (last col) = Infliximab
meta_inflix_CD <- meta_polished[
  toupper(trimws(meta_polished[, 2])) == "CD" &
    toupper(trimws(meta_polished[, ncol(meta_polished)])) == "INFLIXIMAB",
]

# Sample IDs
inflx_CD_ids <- as.character(meta_inflix_CD[, 1])

# Files / names for the CD+Infliximab cohort
Inflix_sampleFiles <- sampleFiles[ sub("_counts.txt$", "", sampleFiles) %in% inflx_CD_ids ]
Inflix_sampleNames <- sub("_counts.txt$", "", Inflix_sampleFiles)

# Conditions (col 3 = Remission / Non-remission) → factor R / NR
status_map <- setNames(as.character(meta_inflix_CD[, 3]), as.character(meta_inflix_CD[, 1]))
Inflix_sampleCondition <- factor(ifelse(status_map[Inflix_sampleNames] == "Remission", "R", "NR"),
                                 levels = c("R","NR"))

## --- DGE for Infliximab + CD cohort 
# Build sample table
sample_table_CD_ifx <- data.frame(
  sampleName = Inflix_sampleNames,
  fileName   = Inflix_sampleFiles,
  condition  = Inflix_sampleCondition,
  stringsAsFactors = FALSE
)


# DESeq2 object and run
ddsHTSeq_CD_ifx <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table_CD_ifx,
                                              design = ~ condition)
colData(ddsHTSeq_CD_ifx)$condition <- factor(colData(ddsHTSeq_CD_ifx)$condition,
                                             levels = c("NR","R"))

dds_CD_ifx <- DESeq(ddsHTSeq_CD_ifx)


# QC: transformed counts, correlation heatmap, PCA
vst_CD_ifx <- varianceStabilizingTransformation(dds_CD_ifx, blind = TRUE)
vst_mat_CD_ifx <- assay(vst_CD_ifx)
vst_cor_CD_ifx <- cor(vst_mat_CD_ifx)

pheatmap(vst_cor_CD_ifx)
plotPCA(vst_CD_ifx, intgroup = "condition")

#rerun the analysis without outlier
if (any(grepl("83257$", Inflix_sampleNames))) {
  keep_ix <- !grepl("83257$", Inflix_sampleNames)
  
  Inflix_sampleNames     <- Inflix_sampleNames[keep_ix]
  Inflix_sampleFiles     <- Inflix_sampleFiles[keep_ix]
  Inflix_sampleCondition <- Inflix_sampleCondition[keep_ix]
  
  sample_table_CD_ifx <- data.frame(
    sampleName = Inflix_sampleNames,
    fileName   = Inflix_sampleFiles,
    condition  = Inflix_sampleCondition,
    stringsAsFactors = FALSE
  )
  
  ddsHTSeq_CD_ifx <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table_CD_ifx,
                                                design = ~ condition)
  colData(ddsHTSeq_CD_ifx)$condition <- factor(colData(ddsHTSeq_CD_ifx)$condition,
                                               levels = c("NR","R"))
  dds_CD_ifx <- DESeq(ddsHTSeq_CD_ifx)
  
  vst_CD_ifx     <- varianceStabilizingTransformation(dds_CD_ifx, blind = TRUE)
  vst_mat_CD_ifx <- assay(vst_CD_ifx)
  vst_cor_CD_ifx <- cor(vst_mat_CD_ifx)
  
  pheatmap(vst_cor_CD_ifx)
  plotPCA(vst_CD_ifx, intgroup = "condition")
}
## PCA seems to be okish, I will proceed with the rest of the analysis.

# Convert Ensembl IDs (human) to HGNC symbols for the CD + Infliximab cohort

library(AnnotationDbi)
library(org.Hs.eg.db)

ensembl_gene_ids_CD_ifx <- sub("\\.\\d+$", "", rownames(dds_CD_ifx))  # remove Ensembl version suffixes
hgnc_symbols_CD_ifx <- mapIds(org.Hs.eg.db,
                              keys      = ensembl_gene_ids_CD_ifx,
                              column    = "SYMBOL",
                              keytype   = "ENSEMBL",
                              multiVals = "first")

mcols(dds_CD_ifx)$hgnc_symbol <- hgnc_symbols_CD_ifx

normalized_counts_CD_ifx <- counts(dds_CD_ifx, normalized = TRUE)
rownames(normalized_counts_CD_ifx) <- ifelse(!is.na(hgnc_symbols_CD_ifx) & hgnc_symbols_CD_ifx != "",
                                             hgnc_symbols_CD_ifx,
                                             rownames(dds_CD_ifx))
View(normalized_counts_CD_ifx)

# Results table (CD + Infliximab) with HGNC symbols for labeling
res_CD_ifx <- results(dds_CD_ifx)
res_CD_ifx$hgnc_symbol <- hgnc_symbols_CD_ifx[rownames(res_CD_ifx)]

res_CD_ifx
summary(res_CD_ifx)

# Filter by adjusted p-value
res05_CD_ifx <- subset(res_CD_ifx, !is.na(padj) & padj < 0.05)
head(res05_CD_ifx)

# Order by padj
resOrdered_CD_ifx <- res05_CD_ifx[order(res05_CD_ifx$padj), ]

# Save (includes symbol column)
write.csv(as.data.frame(resOrdered_CD_ifx), "results-DESeq2-Infliximab_CD.csv", row.names = TRUE)

# MA plot
plotMA(res_CD_ifx, main = "Infliximab_CD_RNA-seq")

# EnhancedVolcano with gene symbols as point labels
library(ggrepel)
library(EnhancedVolcano)

gene_labels_CD_ifx <- ifelse(is.na(res_CD_ifx$hgnc_symbol) | res_CD_ifx$hgnc_symbol == "",
                             rownames(res_CD_ifx), res_CD_ifx$hgnc_symbol)

EnhancedVolcano(as.data.frame(res_CD_ifx),
                lab = gene_labels_CD_ifx,
                x = "log2FoldChange",
                y = "pvalue",
                pCutoff = 1e-3,
                xlim = c(-5, 5))

############
## -------- AutoCore enrichment for Infliximab + CD (SYMBOL-based) --------
library(dplyr)
library(tidyr)
library(clusterProfiler)

# --- Paths ---
autocore_clusters_path <- "/Users/sara/Desktop/Autocore/IBD/Enrichment_Test/Modified_Autocore_Clusters_GeneSymbols copy.txt"
pc_universe_path       <- "/Users/sara/Desktop/Autocore/IBD/Enrichment_Test/symbols2entrezIDs (1).txt"

# 1) DEGs in SYMBOL space (padj < 0.05; adjust if desired)
res_tbl <- as.data.frame(res_CD_ifx) |> tibble::rownames_to_column("ensembl")
stopifnot("hgnc_symbol" %in% names(res_tbl))
res_tbl <- dplyr::rename(res_tbl, symbol = hgnc_symbol)

deg_all_symbols  <- res_tbl |> filter(!is.na(padj), padj < 0.05) |> pull(symbol) |> unique() |> na.omit()
deg_up_symbols   <- res_tbl |> filter(!is.na(padj), padj < 0.05, log2FoldChange > 0) |> pull(symbol) |> unique() |> na.omit()
deg_down_symbols <- res_tbl |> filter(!is.na(padj), padj < 0.05, log2FoldChange < 0) |> pull(symbol) |> unique() |> na.omit()

# 2) Protein-coding universe (SYMBOL) from your file
pc_tbl <- read.delim(pc_universe_path, check.names = TRUE, stringsAsFactors = FALSE)
names(pc_tbl) <- make.names(names(pc_tbl), unique = TRUE)  # normalize headers once

sym_col <- intersect(c("Approved.Symbol","Approved_Symbol","SYMBOL","Symbol","symbol"), names(pc_tbl))[1]
loc_col <- intersect(c("Locus.Type","Locus_Type","locus_type"), names(pc_tbl))[1]
stopifnot(length(sym_col) == 1, length(loc_col) == 1)

pc_universe_symbols <- pc_tbl |>
  filter(!is.na(.data[[sym_col]]),
         .data[[sym_col]] != "",
         grepl("gene with protein product", .data[[loc_col]], ignore.case = TRUE)) |>
  pull(all_of(sym_col)) |>
  unique()

# 3) AutoCore clusters (SYMBOL) → TERM2GENE
clusters_raw <- read.delim(autocore_clusters_path, check.names = FALSE, stringsAsFactors = FALSE)

## 0) Clean header names so none are "" or NA
hdr <- names(clusters_raw)
hdr <- trimws(hdr)
hdr[is.na(hdr) | hdr == ""] <- paste0("X", which(is.na(hdr) | hdr == ""))
names(clusters_raw) <- make.names(hdr, unique = TRUE)

## 1) Lock the expected columns (case-insensitive)
nm <- tolower(names(clusters_raw))
symbol_col  <- names(clusters_raw)[which(nm == "symbol")[1]]
cluster_col <- names(clusters_raw)[which(nm == "cluster_name")[1]]
stopifnot(length(symbol_col) == 1, length(cluster_col) == 1)

## 2) First step: remove NAs in these two columns (base R, not dplyr)
clusters_raw <- clusters_raw[ !is.na(clusters_raw[[symbol_col]]) &
                                !is.na(clusters_raw[[cluster_col]]) , , drop = FALSE ]

## 3) Build TERM2GENE (and keep only protein-coding universe symbols)
TERM2GENE <- dplyr::tibble(
  term = trimws(clusters_raw[[cluster_col]]),
  gene = trimws(clusters_raw[[symbol_col]])
) |>
  dplyr::filter(term != "", gene != "") |>
  dplyr::distinct() |>
  dplyr::filter(gene %in% pc_universe_symbols)

## 4) DEG symbol sets from your CD + Infliximab results (assumes res_CD_ifx has 'hgnc_symbol')
res_tbl <- as.data.frame(res_CD_ifx)
stopifnot("hgnc_symbol" %in% names(res_tbl))

deg_all_symbols  <- unique(na.omit(res_tbl$hgnc_symbol[ !is.na(res_tbl$padj) & res_tbl$padj < 0.05 ]))
deg_up_symbols   <- unique(na.omit(res_tbl$hgnc_symbol[ !is.na(res_tbl$padj) & res_tbl$padj < 0.05 & res_tbl$log2FoldChange > 0 ]))
deg_down_symbols <- unique(na.omit(res_tbl$hgnc_symbol[ !is.na(res_tbl$padj) & res_tbl$padj < 0.05 & res_tbl$log2FoldChange < 0 ]))

## 5) Enrichment (ORA) with SYMBOL IDs and protein-coding background
 

enr_autocore_CD_ifx_all  <- enrich_autocore(deg_all_symbols)
enr_autocore_CD_ifx_up   <- enrich_autocore(deg_up_symbols)
enr_autocore_CD_ifx_down <- enrich_autocore(deg_down_symbols)

# (optional) write results if not NULL:
# if (!is.null(enr_autocore_CD_ifx_all))  write.csv(as.data.frame(enr_autocore_CD_ifx_all),  "AutoCore_enrichment_CD_IFX_ALL.csv",  row.names = FALSE)
# if (!is.null(enr_autocore_CD_ifx_up))   write.csv(as.data.frame(enr_autocore_CD_ifx_up),   "AutoCore_enrichment_CD_IFX_UP.csv",   row.names = FALSE)
# if (!is.null(enr_autocore_CD_ifx_down)) write.csv(as.data.frame(enr_autocore_CD_ifx_down), "AutoCore_enrichment_CD_IFX_DOWN.csv", row.names = FALSE)
