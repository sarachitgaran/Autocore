
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

# Keep only Infliximab-treated samples
meta <- read.csv("../IBD_Metadata.csv", stringsAsFactors = FALSE)
meta_polished<- meta[,c(1,14,25,31,32)]

# subset only Infliximab-treated samples
meta_inflix <- meta_polished[meta_polished[,ncol(meta_polished)] == "Infliximab", ]
dim(meta_inflix) #126 samples

# Vector of Infliximab sample IDs
inflx_ids <- as.character(meta_inflix[, 1])

# New variables prefixed with "Inflix", subset to only those samples
Inflix_sampleFiles <- sampleFiles[sub("_counts.txt$", "",sampleFiles) %in% inflx_ids ]
Inflix_sampleNames <- sub("_counts.txt$", "",Inflix_sampleFiles)

# Subset to Vedolizumab-treated samples and create variables
meta_vedo <- meta_polished[ meta_polished[, ncol(meta_polished)] == "Vedolizumab", ]
vedo_ids <- as.character(meta_vedo[, 1])

Vedo_sampleFiles <- sampleFiles[ sub("_counts.txt$", "", sampleFiles) %in% vedo_ids ]
Vedo_sampleNames <- sub("_counts.txt$", "", Vedo_sampleFiles)


## DGE for Infliximab ##

# From meta_inflix (col1 = sample ID, col3 = Remission/Non-remission) to factor R/NR
status_map <- setNames(as.character(meta_inflix[,3]), as.character(meta_inflix[,1]))
Inflix_sampleCondition <- factor(ifelse(status_map[Inflix_sampleNames] == "Remission", "R", "NR"),
                          levels = c("R","NR"))

table <- data.frame(sampleName = Inflix_sampleNames,
                    fileName = Inflix_sampleFiles,
                    condition = Inflix_sampleCondition)

# Creating DESeq2 object
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = table,
                                       design = ~condition)

ddsHTSeq$condition

# Setting the factor levels
treatments <- c("NR","R")

colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,
                                      levels = treatments)

# Run DESeq2 differential expression analysis
dds <- DESeq(ddsHTSeq)

View(counts(dds, normalized=T))

#### Quality Control ####

# Transform counts for data visualization
vst <- varianceStabilizingTransformation(dds, blind=T)

# Extract the rlog matrix from the object and compute pairwise correlation values
vst_mat <- assay(vst)
vst_cor <- cor(vst_mat)


# Plot heatmap (SDM)
library(pheatmap)
pheatmap(vst_cor)

# Plot PCA 
plotPCA(vst, intgroup="condition")

## Seems like there is an outlier- 
## I will remove this sample and redo the analysis so far.

drop_idx <- grepl("83257$", Inflix_sampleNames)

Inflix_sampleNames      <- Inflix_sampleNames[!drop_idx]
Inflix_sampleFiles      <- Inflix_sampleFiles[!drop_idx]
Inflix_sampleCondition  <- Inflix_sampleCondition[!drop_idx]

table <- data.frame(
  sampleName = Inflix_sampleNames,
  fileName   = Inflix_sampleFiles,
  condition  = Inflix_sampleCondition,
  stringsAsFactors = FALSE
)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = table, design = ~ condition)
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition, levels = c("NR","R"))
dds <- DESeq(ddsHTSeq)

# Transform counts for data visualization
vst <- varianceStabilizingTransformation(dds, blind=T)

# Extract the rlog matrix from the object and compute pairwise correlation values
vst_mat <- assay(vst)
vst_cor <- cor(vst_mat)


# Plot heatmap (SDM)
library(pheatmap)
pheatmap(vst_cor)

# Plot PCA 
plotPCA(vst, intgroup="condition")


#convert to gene symbol
# Convert Ensembl IDs (human) to HGNC symbols with informative variable names
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
library(AnnotationDbi)
library(org.Hs.eg.db)

ensembl_gene_ids <- sub("\\.\\d+$", "", rownames(dds))  # remove Ensembl version suffixes
hgnc_symbols <- mapIds(org.Hs.eg.db,
                       keys      = ensembl_gene_ids,
                       column    = "SYMBOL",
                       keytype   = "ENSEMBL",
                       multiVals = "first")

mcols(dds)$hgnc_symbol <- hgnc_symbols

normalized_counts_matrix <- counts(dds, normalized = TRUE)
rownames(normalized_counts_matrix) <- ifelse(!is.na(hgnc_symbols) & hgnc_symbols != "",
                                             hgnc_symbols,
                                             rownames(dds))
View(normalized_counts_matrix)


# Results table will be generated using results() which will include:
# base mean, Log2 fold changes, p values and adjusted p values
# DESeq2 results with HGNC symbols for labeling
res <- results(dds)
res$hgnc_symbol <- hgnc_symbols[rownames(res)]

res
summary(res)

# Filter by adjusted p-value
res05 <- subset(res, !is.na(padj) & padj < 0.05)
head(res05)

# Order by padj
resOrdered <- res05[order(res05$padj),]

# Save (includes symbol column)
write.csv(as.data.frame(resOrdered), "results-DESeq2-normalized.csv")

# MA plot
plotMA(res, main = "RNAseq experiment")

# EnhancedVolcano with gene symbols as point labels
library(ggrepel)
library(EnhancedVolcano)

gene_labels <- ifelse(is.na(res$hgnc_symbol) | res$hgnc_symbol == "",
                      rownames(res), res$hgnc_symbol)

EnhancedVolcano(as.data.frame(res),
                lab = gene_labels,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 10e-4,
                xlim = c(-5, 5))








